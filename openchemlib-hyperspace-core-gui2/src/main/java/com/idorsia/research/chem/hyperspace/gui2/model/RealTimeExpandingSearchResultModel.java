package com.idorsia.research.chem.hyperspace.gui2.model;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.coords.InventorTemplate;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;

public class RealTimeExpandingSearchResultModel {

    private CombinatorialSearchResultModel resultModel;

    private int maxExpandedHits = 32000;

    private ConcurrentHashMap<String, SynthonSpace.CombinatorialHit> assembledMolecules;

    // this one is also used for synchronization of computation
    private ConcurrentHashMap<SynthonSpace.CombinatorialHit, List<String>> assembledMolecules2;

    // this defines the row order of the table model
    // !! ONLY ACCESS / MANIPULATE FROM WITHIN swing event handling threaad !!
    private List<String> moleculeOrder = new ArrayList<>();

    public RealTimeExpandingSearchResultModel(CombinatorialSearchResultModel resultModel, int maxExpandedHits) {
        this.resultModel = resultModel;
        this.maxExpandedHits = maxExpandedHits;

        this.assembledMolecules  = new ConcurrentHashMap<>();
        this.assembledMolecules2 = new ConcurrentHashMap<>();

        this.initExpansionThreadPool();
        this.initCoordinateInventor();

        this.resultModel.addListener(new CombinatorialSearchResultModel.CombinatorialSearchResultModelListener() {
            @Override
            public void resultsChanged() {
                processResultsChanged();
            }
        });
        this.processResultsChanged();
    }

    private RealTimeExpandingSearchResultModel getThis() {
        return this;
    }

    private final List<Future> expansionTasks = Collections.synchronizedList(new ArrayList<>());

    private void processResultsChanged() {
        // we may end up here from edt.. therefore do this
        // async..
        Thread tri = new Thread() {
            @Override
            public void run() {
                synchronized(assembledMolecules2) {
                    long tsa = System.currentTimeMillis();
                    //System.out.println("monitor asm2 entered");
                    List<SynthonSpace.CombinatorialHit> all_hits = new ArrayList<>(resultModel.getHits());
                    all_hits.removeAll( assembledMolecules2.keySet() );
                    for(SynthonSpace.CombinatorialHit chi : all_hits) {
                        // immediately put placeholder for molecules, such that the line above works..
                        assembledMolecules2.put(chi,new ArrayList<>());
                        SynthonSpace.CombinatorialHit fchi = chi;
                        Runnable ri = new Runnable() {
                            @Override
                            public void run() {
                                if(moleculeOrder.size()>maxExpandedHits) {return;}
                                List<SynthonAssembler.ExpandedCombinatorialHit> exp_hits = SynthonAssembler.expandCombinatorialHit(fchi,1024);
                                List<String> new_molecules = new ArrayList<>();
                                synchronized(assembledMolecules2) {
                                    if(moleculeOrder.size()>maxExpandedHits) {return;}
                                    List<String> processed = new ArrayList<>();
                                    for(SynthonAssembler.ExpandedCombinatorialHit xi : exp_hits) {
                                        if(Thread.currentThread().isInterrupted()) {break;}
                                        String processed_idcode = processResultStructure(xi.assembled_idcode);
                                        assembledMolecules.put(processed_idcode,fchi);
                                        processed.add(processed_idcode);
                                    }
                                    //new_molecules = new ArrayList<>(exp_hits.stream().map(ci -> ci.assembled_idcode).collect(Collectors.toList()));
                                    new_molecules = processed;
                                    assembledMolecules2.put(fchi,
                                            new_molecules);
                                }
                                List<String> final_new_molecules = new_molecules;
                                SwingUtilities.invokeLater(new Runnable(){
                                    @Override
                                    public void run() {
                                        //System.out.println("SwingEDT: enter");
                                        long tsa = System.currentTimeMillis();
                                        if(moleculeOrder.size()>maxExpandedHits) {return;}
                                        // update internal table model..
                                        int size_old = moleculeOrder.size();
                                        moleculeOrder.addAll( new ArrayList<>( final_new_molecules ) );
                                        int size_new = moleculeOrder.size();
                                        tableModel.fireTableRowsInserted(size_old,size_new);
                                        //System.out.println("SwingEDT: leave, time= "+(System.currentTimeMillis()-tsa));
                                    }
                                });
                                fireResultsChanged();
                            }
                        };

                        if(moleculeOrder.size()<maxExpandedHits) {
                            if(getThis().expansionThreadPool==null || getThis().expansionThreadPool.isShutdown() ) {
                                getThis().restartThreadpool();
                            }
                            Future fExp = expansionThreadPool.submit(ri);
                            synchronized(expansionTasks) {
                                expansionTasks.add(fExp);
                            }
                        }
                        //ri.setPriority(Thread.MIN_PRIORITY);
                        //ri.start();
                    }
                    //System.out.println("monitor asm2 left, time= "+(tsa-System.currentTimeMillis()));
                }

            }
        };
        tri.start();
    }

    // Create a ThreadPoolExecutor with 4 threads
    private ThreadPoolExecutor expansionThreadPool;

    private void initExpansionThreadPool() {
        if(!expansionTasks.isEmpty()) {
            List<Future> allTasks = new ArrayList<>();
            synchronized(expansionTasks) {
                allTasks = new ArrayList<>(expansionTasks);
            }
            for(Future fi : allTasks) {
                fi.cancel(true);
                expansionTasks.remove(fi);
            }
        }

        expansionThreadPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(2);

        // Set the lowest priority for each thread
        expansionThreadPool.setThreadFactory(new ThreadFactory() {
            @Override
            public Thread newThread(Runnable r) {
                Thread t = new Thread(r);
                t.setPriority(Thread.MIN_PRIORITY);
                return t;
            }
        });
    }

    private boolean highlightSubstructure = true;
    private boolean alignSubstructure     = false;

    public void setStructurePostprocessOptions( boolean highlightSubstructure, boolean alignSubstructure ) {
        boolean reprocessNeeded = false;
        if(this.alignSubstructure != alignSubstructure) {
            this.alignSubstructure = alignSubstructure;
            initCoordinateInventor();
            reprocessNeeded = true;
        }
        if(this.highlightSubstructure != highlightSubstructure) {
            this.highlightSubstructure = highlightSubstructure;
            reprocessNeeded = true;
        }
        if(reprocessNeeded) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    moleculeOrder.clear();
                    assembledMolecules.clear();
                    assembledMolecules2.clear();
                }
            });
            processResultsChanged();
        }
    }

    public boolean isHighlightSubstructure() {
        return highlightSubstructure;
    }

    public boolean isAlignSubstructure() {
        return alignSubstructure;
    }

    private CoordinateInventor coordinateInventor;
    private void initCoordinateInventor() {
        DescriptorHandlerLongFFP512 ffp = new DescriptorHandlerLongFFP512();
        coordinateInventor = new CoordinateInventor();
        if( this.resultModel.getQuery()!=null && alignSubstructure ) {
            StereoMolecule qi = this.resultModel.getQuery();
            InventorTemplate it = new InventorTemplate(qi, ffp.createDescriptor(qi), true);
            coordinateInventor.setCustomTemplateList(Collections.singletonList(it));
        }
    }

    public CombinatorialSearchResultModel getCombinatorialSearchResultModel() {
        return this.resultModel;
    }

    private List<RealTimeExpandingSearchResultModelListener> listeners = new ArrayList<>();

    private RealTimeExpandingTableModel tableModel = new RealTimeExpandingTableModel();
    public RealTimeExpandingTableModel getTableModel() {
        return tableModel;
    }

    public void addListener(RealTimeExpandingSearchResultModelListener li) {
        this.listeners.add(li);
    }

    public boolean removeListener(RealTimeExpandingSearchResultModelListener li) {
        return this.listeners.remove(li);
    }

    private void fireResultsChanged() {
        for(RealTimeExpandingSearchResultModelListener li : listeners) {
            li.resultsChanged();
        }
    }


    public String processResultStructure(String mol) {

        StereoMolecule mi = HyperspaceUtils.parseIDCode(mol);
        coordinateInventor.invent(mi);
        if(this.highlightSubstructure) {
            HyperspaceUtils.setHighlightedSubstructure(mi, resultModel.getQuery());
        }
        Canonizer ci = new Canonizer(mi);
        String result = ci.getIDCode() +" "+ci.getEncodedCoordinates();
        return result;
    }

    public String getResultsInfoString() {
        long resultsTotal = resultModel.getHits().stream().mapToLong( xi ->
                xi.hit_fragments.values().stream().mapToLong(x->(long)x.size()).reduce((x,y)->x*y).getAsLong()).sum();

        long numShown = 0;
        //synchronized(assembledMolecules2) {
        numShown = this.moleculeOrder.size();
        //}
        String ri = String.format("Results: %6d  Showing: %6d", resultsTotal, numShown);
        return ri;
    }


    public static interface RealTimeExpandingSearchResultModelListener {
        public void resultsChanged();
    }

    public class RealTimeExpandingTableModel extends AbstractTableModel {

        public String getStructureData(int row) {
            return moleculeOrder.get(row);
        }

        @Override
        public String getColumnName(int column) {
            switch(column){
                case 0: return "Structure";
            }
            return null;
        }

        @Override
        public int getRowCount() {
            return moleculeOrder.size();
        }

        @Override
        public int getColumnCount() {
            return 1;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            return moleculeOrder.get(rowIndex);
        }
    }

    public void shutdownThreadpool() {
        if(!expansionTasks.isEmpty()) {
            List<Future> allTasks = new ArrayList<>(expansionTasks);
            for(Future fi : allTasks) {
                fi.cancel(true);
                expansionTasks.remove(fi);
            }
        }
        this.expansionThreadPool.shutdown();
    }

    public void restartThreadpool() {
        this.expansionThreadPool.shutdown();
        this.initExpansionThreadPool();
    }

    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        shutdownThreadpool();
    }
}
