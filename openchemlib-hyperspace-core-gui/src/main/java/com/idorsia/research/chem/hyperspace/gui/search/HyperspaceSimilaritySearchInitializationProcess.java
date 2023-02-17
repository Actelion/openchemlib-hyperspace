package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpace3;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import org.apache.commons.io.input.CountingInputStream;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

public class HyperspaceSimilaritySearchInitializationProcess extends AbstractHyperspaceProcess implements AbstractHyperspaceProcess.HasProgress {

    private long file_size_in_bytes;
    private CountingInputStream counting_in_stream = null;

    private String file;

    private double sim_search_initialized = 0.0;
    private double ratio_complete = 0.0;


    HyperspaceSimilaritySearch search_provider;
    SynthonSimilaritySpace3 space;


    public HyperspaceSimilaritySearchInitializationProcess(AbstractSearchProvider search_provider, String file) {
        this.search_provider = (HyperspaceSimilaritySearch) search_provider;
        this.file = file;
        setProcessStatus(ProcessStatus.WAITING);
    }

    public void startInitializationAsync() {
        setProcessStatus(ProcessStatus.COMPUTING);
        Thread t_init = new Thread() {
            @Override
            public void run() {

                try {
                    File fi_file = new File(file);
                    file_size_in_bytes = fi_file.length();

                    FileInputStream file_in = new FileInputStream(fi_file);
                    //counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
                    counting_in_stream = new CountingInputStream( new BufferedInputStream(file_in));

                    ObjectInputStream in = new ObjectInputStream( new GZIPInputStream(counting_in_stream) );


                    //SynthonSpace space_a = null;
                    //space_a = (SynthonSpace) in.readObject();
                    //space_a.initAfterJavaDeserialization();

                    SynthonSimilaritySpace3 space_a = null;
                    space_a = (SynthonSimilaritySpace3) in.readObject();
                    space_a.initAfterJavaDeserialization();

                    //space = space_a;

                    System.out.println("Loaded space: " + space_a.getSpaceInfoString());

                    if(false) {
                        System.out.println("Now initialize synthon topology data");

                        Consumer<Double> init_progress = new Consumer<Double>() {
                            @Override
                            public void accept(Double aDouble) {
                                sim_search_initialized = aDouble;
                            }
                        };

                        SynthonSimilaritySpace3 space_3 = new SynthonSimilaritySpace3(space_a);
                        space_3.initFastSimilaritySearchers(6, init_progress);

                        space = space_3;
                    }
                    else {
                        space = space_a;
                    }

                    //setStatus(AbstractSearchProvider.SearchProviderStatus.READY); // we do this in the SearchProvider..
                    search_provider.setSpace(space);
                    setProcessStatus(ProcessStatus.DONE);

                    return;

                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (ClassNotFoundException e) {
                    e.printStackTrace();
                }
                setProcessStatus(ProcessStatus.FAILED);
                //setStatus(AbstractSearchProvider.SearchProviderStatus.ERROR); // we do this in the SearchProvider..
            }
        };
        setProcessStatus(ProcessStatus.COMPUTING);
        t_init.start();
    }

    public double getStatusOfInitialization() {
        //if(status== AbstractSearchProvider.SearchProviderStatus.READY) {return 1.0;}
        if (this.getProcessStatus() == ProcessStatus.DONE) {
            return 1.0;
        }
        if (counting_in_stream == null) {
            return 0.0;
        }
        // OLD: all read means 30 percent done, when all initialized then 90%
        // all read means 95% done.
        double ratio_done = 0.95 * counting_in_stream.getByteCount() / file_size_in_bytes;
        //double ratio_done_2 = 0.6* this.sim_search_initialized;

        return ratio_done;
    }

    @Override
    public String getName() {
        return "Load Synthon Space";
    }


    public SynthonSpace getSpace() {
        return this.space;
    }

    //@Override
    public List<StereoMolecule> getQueryStructures() {
        return new ArrayList<>();
    }

    @Override
    public double getProgress() {
        return this.getStatusOfInitialization();
    }
}