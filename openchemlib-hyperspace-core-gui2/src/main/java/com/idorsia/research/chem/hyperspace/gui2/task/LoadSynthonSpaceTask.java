package com.idorsia.research.chem.hyperspace.gui2.task;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;
import com.idorsia.research.chem.hyperspace.gui2.model.LoadedSynthonSpace;
import org.apache.commons.io.input.CountingInputStream;

import javax.swing.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

public class LoadSynthonSpaceTask extends SwingWorker<SynthonSpace,Double> implements HyperspaceTask {

    private LeetHyperspaceModel model;
    //private LeetHyperspaceView  view;

    private String filepath;

    public LoadSynthonSpaceTask(LeetHyperspaceModel model, String filepath) {
        this.model = model;
        //this.view = view;
        this.filepath = filepath;
    }

    @Override
    public String getName() {
        return "Load space "+ (new File(filepath)).getName();
    }

    public SwingWorker getThisWorker() {return this;}

    @Override
    protected SynthonSpace doInBackground() throws Exception {
        SynthonSpace space = null;
        //CountingInputStream counting_in_stream = null;


        Thread ti = null;

        try {
            File fi_file = new File(filepath);
            long file_size_in_bytes = fi_file.length();

            FileInputStream file_in = new FileInputStream(fi_file);
            //counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
            //counting_in_stream = new CountingInputStream( new GZIPInputStream( new BufferedInputStream(file_in)));
            CountingInputStream counting_in_stream = new CountingInputStream( new BufferedInputStream(file_in));
            ObjectInputStream in = new ObjectInputStream( new GZIPInputStream( counting_in_stream ) );


            ti = new Thread() {
                @Override
                public void run() {
                    while( !this.isInterrupted() ) {
                        double ratio_done = 0.0;
                        if(getThisWorker().isDone()) {ratio_done = 1.0;}
                        if(counting_in_stream==null) {
                            ratio_done = 0.0;
                        }
                        // all read means 80 percent done:y
                        ratio_done = 0.8 * counting_in_stream.getByteCount() / file_size_in_bytes;
                        publish( ratio_done );
                    }
                }
            };
            ti.start();


            SynthonSpace space_a = null;
            space_a = (SynthonSpace) in.readObject();
            space_a.initAfterJavaDeserialization();
            space = space_a;

            System.out.println("Loaded space: "+space.getSpaceInfoString());

            if(true) {
                for( String rxnid : space.getRxnIds()) {
                    SynthonSpace finalSpace = space;
                    List<SynthonSpace.FragType> ssets = new ArrayList<>(space.getFragTypes(rxnid).values());
                    long size = ssets.stream().mapToLong(xi -> finalSpace.getSynthonSet(rxnid,xi.frag).size() ).reduce( (x,y) -> x*y ).getAsLong();
                    boolean billionClubRxn = (size >= 1e9);
                    System.out.println(rxnid+" -> " + ((billionClubRxn)?"[!!BILLION+!!]":"") + size + " : "+ssets.stream().map(xi -> ""+finalSpace.getSynthonSet(rxnid,xi.frag).size() ).collect(Collectors.joining("x")));
                }
            }

            //setStatus(AbstractSearchProvider.SearchProviderStatus.READY); // we do this in the SearchProvider..
            //setProcessStatus(AbstractHyperspaceProcess.ProcessStatus.DONE);

            //return space;

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(ti!=null) {
                ti.interrupt();
            }
        }
        return space;
    }

    @Override
    protected void process(List<Double> chunks) {
        // Update the GUI with the progress values received
        double latestProgress = chunks.get(chunks.size() - 1);
        this.setProgress(  Math.min( 100, Math.max(0, (int) (100.0*latestProgress))) );
    }

    @Override
    protected void done() {
        try {
            this.setProgress(100);
            String name = (new File(this.filepath)).getName();
            model.addSynthonSpace( new LoadedSynthonSpace(this.get(),name));
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        } catch (ExecutionException e) {
            throw new RuntimeException(e);
        }
    }
}
