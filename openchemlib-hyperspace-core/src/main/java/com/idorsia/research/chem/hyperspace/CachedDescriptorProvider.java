package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerStandard2DFactory;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class CachedDescriptorProvider {

    //SynthonSpace space = null;
    String mDescriptorHandlerShortName;

    public CachedDescriptorProvider( String descriptorhandler_shortName ) {
        this.mDescriptorHandlerShortName = descriptorhandler_shortName;
    }

    public ConcurrentHashMap<String, BitSet> cache_FPs = new ConcurrentHashMap<>();

    /**
     * Now Thread-safe!
     *
     * @param mi
     * @return
     */
    public BitSet getFP_cached(StereoMolecule mi) {
        String m_idcode = null;
        synchronized(mi) {
            m_idcode = mi.getIDCode();
        }
        if(cache_FPs.containsKey(m_idcode)) {
            return (BitSet) cache_FPs.get(m_idcode).clone();
        }
        BitSet fpi = SynthonSpace.getFP( this.getDescriptorHandler_PerThread() , mi );
        cache_FPs.put(m_idcode,fpi);
        return (BitSet) fpi.clone();
    }



    DescriptorHandler<long[],StereoMolecule> mFP = null;
    transient Map<Thread,DescriptorHandler<long[],StereoMolecule>> mFPs_forThreads = new HashMap<>();

    public synchronized DescriptorHandler<long[],StereoMolecule> getDescriptorHandler_PerThread() {
        if(mFP==null) {
            // parse and restore the descriptor handler (based on mDescriptorHandlerShortName)
            this.mFP = resolveDescriptorHandlerFromName(mDescriptorHandlerShortName);
        }

        DescriptorHandler<long[],StereoMolecule> dhi = mFPs_forThreads.get( Thread.currentThread() );
        if(dhi==null) {
            dhi = mFP.getThreadSafeCopy();
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
