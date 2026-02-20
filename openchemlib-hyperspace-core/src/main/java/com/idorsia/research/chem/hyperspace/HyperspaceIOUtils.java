package com.idorsia.research.chem.hyperspace;

import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import org.apache.commons.io.input.CountingInputStream;

import java.io.*;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class HyperspaceIOUtils {

    public static SynthonSpace loadSynthonSpace(String file) throws IOException {
        //try {
        File fi_file = new File(file);
        long file_size_in_bytes = fi_file.length();

        FileInputStream file_in = new FileInputStream(fi_file);
        CountingInputStream counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
        //ObjectInputStream in = new ObjectInputStream(counting_in_stream);
        ObjectInputStream in = new ObjectInputStream( new GZIPInputStream( counting_in_stream ) );
        SynthonSpace space = null;

        Thread t_progress_monitor = new Thread() {
            @Override
            public void run() {
                while (true) {
                    try {
                        Thread.sleep(5000);
                    } catch (InterruptedException e) {
                        System.out.println("monitor done!");
                        break;
                    }
                    long bytes_read = counting_in_stream.getByteCount();
                    double ratio = (1.0 * bytes_read) / file_size_in_bytes;
                    System.out.println("read: " + String.format("%4.2f", (ratio * 100)) + " i.e. " + (bytes_read / 1024.0) + " kb / " + (file_size_in_bytes / 1024.0) + " kb");
                }
            }
        };
        t_progress_monitor.start();

        try {
            space = (SynthonSpace) in.readObject();
        } catch (OptionalDataException ex) {
            ex.printStackTrace();
        }
        catch (ClassNotFoundException ex) {
            ex.printStackTrace();
        }

        t_progress_monitor.interrupt();
        space.initAfterJavaDeserialization();

        return space;
    }

    public static void saveSynthonSpace(SynthonSpace sspace, String file) {
        try {
            System.out.println("Start writing the output!");
            File f = new File(file);
            System.out.println("Into File: "+f.getName());
            if(true) {
                //ObjectOutputStream out = new ObjectOutputStream( new BufferedOutputStream( new FileOutputStream(f) ) );
                ObjectOutputStream out = new ObjectOutputStream( new GZIPOutputStream( new BufferedOutputStream( new FileOutputStream(f) ) ) );
                out.writeObject(sspace);
                out.flush();
                out.close();
                System.out.println("Done!");
            }
        }
        catch(IOException ex) {
            ex.printStackTrace();
        }
    }

    public static DownsampledSynthonSpace loadDownsampledSynthonSpace(String file) throws IOException {
        try (ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(file))))) {
            return (DownsampledSynthonSpace) in.readObject();
        } catch (ClassNotFoundException e) {
            throw new IOException("Failed to deserialize DownsampledSynthonSpace", e);
        }
    }

    public static RawSynthonSpace loadRawSynthonSpace(String file) throws IOException {
        return RawSynthonSpaceIO.read(Path.of(file));
    }

    public static void saveRawSynthonSpace(RawSynthonSpace space, String file) throws IOException {
        RawSynthonSpaceIO.write(space, Path.of(file));
    }

}
