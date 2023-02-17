package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import org.apache.commons.io.input.CountingInputStream;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class HyperspaceSynthonSpaceInitializationProcess extends AbstractHyperspaceProcess implements AbstractHyperspaceProcess.HasProgress {

    private long file_size_in_bytes;
    private CountingInputStream counting_in_stream = null;

    private String file;
    private double ratio_complete = 0.0;

    AbstractSearchProvider search_provider;
    SynthonSpace           space;


    public HyperspaceSynthonSpaceInitializationProcess(AbstractSearchProvider search_provider, String file) {
        this.search_provider = search_provider;
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
                    //counting_in_stream = new CountingInputStream( new GZIPInputStream( new BufferedInputStream(file_in)));
                    counting_in_stream = new CountingInputStream( new BufferedInputStream(file_in));

                    ObjectInputStream in = new ObjectInputStream( new GZIPInputStream( counting_in_stream ) );

                    SynthonSpace space_a = null;
                    space_a = (SynthonSpace) in.readObject();
                    space_a.initAfterJavaDeserialization();
                    space = space_a;

                    System.out.println("Loaded space: "+space.getSpaceInfoString());

                    //setStatus(AbstractSearchProvider.SearchProviderStatus.READY); // we do this in the SearchProvider..
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
        if(this.getProcessStatus()==ProcessStatus.DONE) {return 1.0;}
        if(counting_in_stream==null) {
            return 0.0;
        }
        // all read means 80 percent done:y
        double ratio_done = 0.8 * counting_in_stream.getByteCount() / file_size_in_bytes;
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