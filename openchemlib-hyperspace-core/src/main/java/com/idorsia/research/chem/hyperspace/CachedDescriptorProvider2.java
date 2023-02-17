package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerStandard2DFactory;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Differences to CachedDescriptorProvider:
 *
 * When you get the Object, then you get the actual
 * reference to the stored object, and not a new equivalent
 * object.
 *
 *
 *
 */
public class CachedDescriptorProvider2<T> {

    //SynthonSpace space = null;
    String mDescriptorHandlerShortName;

    public CachedDescriptorProvider2( String descriptorhandler_shortName ) {
        this.mDescriptorHandlerShortName = descriptorhandler_shortName;
    }


    public ConcurrentHashMap<String, T> cache_FPs = new ConcurrentHashMap<>();

    /**
     * Now Thread-safe!
     *
     * @param mi
     * @return
     */
    public T getFP_cached(StereoMolecule mi) {
        String m_idcode = null;
        synchronized(mi) {
            m_idcode = mi.getIDCode();
        }
        if(cache_FPs.containsKey(m_idcode)) {
            return (T) cache_FPs.get(m_idcode);
        }
        T fpi = this.getDescriptorHandler_PerThread().createDescriptor(mi);
        cache_FPs.put(m_idcode,fpi);
        return (T) fpi;
    }



    DescriptorHandler<long[],StereoMolecule> mFP = null;
    transient Map<Thread,DescriptorHandler<T,StereoMolecule>> mFPs_forThreads = new HashMap<>();

    public synchronized DescriptorHandler<T,StereoMolecule> getDescriptorHandler_PerThread() {
        if(mFP==null) {
            // parse and restore the descriptor handler (based on mDescriptorHandlerShortName)
            this.mFP = resolveDescriptorHandlerFromName(mDescriptorHandlerShortName);
        }

        DescriptorHandler<T,StereoMolecule> dhi = mFPs_forThreads.get( Thread.currentThread() );
        if(dhi==null) {
            dhi = (DescriptorHandler<T,StereoMolecule>) mFP.getThreadSafeCopy();
            mFPs_forThreads.put(Thread.currentThread(),dhi);
        }

        return dhi;
    }

    private static DescriptorHandler<long[],StereoMolecule> resolveDescriptorHandlerFromName(String shortName) {
        if(shortName.equals("FFP1024_plus_ffp")) {
            return new DescriptorHandlerLongFFP1024_plus("ffp");
        }
        else if(shortName.equals("FFP1024_plus_pfp")) {
            return new DescriptorHandlerLongFFP1024_plus("pfp");
        }
        else if(shortName.equals("PPC2048")) {
            return new DescriptorHandlerPPCore();
        }
        else {
            return DescriptorHandlerStandard2DFactory.getFactory().create(shortName);
        }
    }


}
