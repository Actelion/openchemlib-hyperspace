package com.idorsia.research.chem.hyperspace.gui.process;

import com.actelion.research.gui.VerticalFlowLayout;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import info.clearthought.layout.TableLayout;
import info.clearthought.layout.TableLayoutConstants;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.ListDataEvent;
import javax.swing.event.ListDataListener;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.List;

public class JProcessList extends JPanel {

    public static class ProcessListModel extends AbstractListModel<AbstractHyperspaceProcess> {

        private List<AbstractHyperspaceProcess> processes = new ArrayList<>();

        public void addProcess(AbstractHyperspaceProcess p) {
            this.processes.add(p);
            this.fireContentsChanged(this,this.processes.size()-1,this.processes.size());
        }
        public void removeProcess(AbstractHyperspaceProcess p) {
            this.processes.remove(p);
            this.fireContentsChanged(this,this.processes.size()-1,this.processes.size());
        }

        public List<AbstractHyperspaceProcess> getProcessList() {
            return new ArrayList<>(this.processes);
        }

        @Override
        public int getSize() {
            return processes.size();
        }

        @Override
        public AbstractHyperspaceProcess getElementAt(int index) {
            return processes.get(index);
        }
    }

    ProcessListModel model;



    // GUI
    JScrollPane scp_a;
    JPanel      jp_a;

    boolean autoScrollDown = true;

    // GUI end

    public JProcessList(ProcessListModel model) {
        this.model = model;
        this.initGUI();

        this.model.addListDataListener(new ListDataListener() {
            @Override
            public void intervalAdded(ListDataEvent e) {
                updateList();
            }
            @Override
            public void intervalRemoved(ListDataEvent e) {
                updateList();
            }
            @Override
            public void contentsChanged(ListDataEvent e) {
                updateList();
            }
        });

        // init update thread:
        Thread t_update = new Thread() {
            @Override
            public void run() {
                while(true) {
                    try {
                        Thread.sleep(250);
                        SwingUtilities.invokeLater(new Runnable() {
                            @Override
                            public void run() {
                                for( Component ci : jp_a.getComponents() ) {
                                    if(ci instanceof JProcessPanel) {
                                        ((JProcessPanel)ci).updatePanelGUI();
                                    }
                                }
                            }
                        });

                    } catch (InterruptedException e) {
                        System.out.println("interrupted -> exit monitor thread");
                        break;
                    }
                }
            }
        };
        t_update.start();
    }

    public void initGUI() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.jp_a  = new JPanel();
        this.scp_a = new JScrollPane(this.jp_a);

        this.add(this.scp_a,BorderLayout.CENTER);

        this.setBorder( BorderFactory.createTitledBorder("Processes") );
    }

    @Override
    public Dimension getMaximumSize() {
        return new Dimension( HiDPIHelper.scale(140) , 10000 );
    }

    public void updateList() {
        this.jp_a.removeAll();
        this.jp_a.setLayout(new VerticalFlowLayout(VerticalFlowLayout.LEFT,VerticalFlowLayout.TOP));
        for(AbstractHyperspaceProcess pi : this.model.getProcessList()) {
            this.jp_a.add(new JProcessPanel(pi));
        }

        if(autoScrollDown) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    JScrollBar vertical = scp_a.getVerticalScrollBar();
                    vertical.setValue( vertical.getMaximum() );
                }
            });
        }
    }

    public class JProcessPanel extends JPanel {
        final AbstractHyperspaceProcess p;
        JPanel jp;
        JLabel jl_a;
        JTextField jtf_b;
        JTextField jtf_c;
        JProgressBar jpb_a = null;
        public JProcessPanel(AbstractHyperspaceProcess p) {
            this.p = p;
            initPanelGUI();

            //Object height_constraint = HiDPIHelper.scale(20);
            //Object height_constraint = TableLayout.PREFERRED;
            //this.setMaximumSize(new Dimension( HiDPIHelper.scale(200) , height_constraint ));
            //this.setPreferredSize(new Dimension( HiDPIHelper.scale(200) , HiDPIHelper.scale(20) ));
            //this.setBorder(new LineBorder(Color.red));
            this.setBorder(new EtchedBorder());
        }
        private void initPanelGUI() {
            this.removeAll();
            //this.setLayout(new BorderLayout());
            //this.setLayout(new TableLayout(new double[][]{ { 0.75,0.25} , {HiDPIHelper.scale(120)} }));
            //this.setLayout(new TableLayout(new double[][]{ { TableLayoutConstants.FILL,HiDPIHelper.scale(60)} , {HiDPIHelper.scale(40)} }));
            this.setLayout(new TableLayout(new double[][]{ { TableLayoutConstants.FILL,HiDPIHelper.scale(120)} , {HiDPIHelper.scale(20)} }));
            this.jp = new JPanel();
            //this.jp.setLayout(new FlowLayout(FlowLayout.LEFT));
            TableLayout tb_a = new TableLayout(new double[][]{ { TableLayoutConstants.PREFERRED , TableLayoutConstants.FILL , TableLayoutConstants.PREFERRED } , {HiDPIHelper.scale(20)} });
            tb_a.setHGap(2);
            this.jp.setLayout(tb_a);
            //this.jp.setLayout(new TableLayout(new double[][]{ {0.85,0.15} , {1.0} }));
            this.jl_a = new JLabel("?");

            if(p instanceof AbstractHyperspaceSearchProcess) {
                AbstractHyperspaceSearchProcess hsp = (AbstractHyperspaceSearchProcess) p;
                JPopupMenu jmi = new JPopupMenu();
                Action presentResults = new AbstractAction("Show Results") {
                    @Override
                    public void actionPerformed(ActionEvent actionEvent) {
                        if(hsp.getProcessStatus()== AbstractHyperspaceProcess.ProcessStatus.DONE) {
                            (hsp).getSearchProvider().presentSearchResults(hsp.getSearchResults());
                        }
                    }
                };
                hsp.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
                    @Override
                    public void processStatusChanged() {
                        if(hsp.getProcessStatus()== AbstractHyperspaceProcess.ProcessStatus.DONE) {
                            presentResults.setEnabled(true);
                            jl_a.setToolTipText(hsp.getSearchConfiguration().getQueryMolecules().get(0).getIDCode()+" Results:"+hsp.getSearchResults().size()+":"+SubstructureSearchHelper.countHits(hsp.getSearchResults()) );
                        }
                    }
                });
                presentResults.setEnabled( hsp.getProcessStatus()==AbstractHyperspaceProcess.ProcessStatus.DONE );
                jmi.add(new JMenuItem(presentResults));
                this.jl_a.setComponentPopupMenu(jmi);

                try {
                    this.jl_a.setToolTipText(hsp.getSearchConfiguration().getQueryMolecules().get(0).getIDCode()+ " Results:"+((hsp.getSearchResults()!=null)?""+hsp.getSearchResults().size()+":"+":"+SubstructureSearchHelper.countHits(hsp.getSearchResults()):"?"));
                }
                catch(Exception ex) {
                }

            }

            this.jtf_b = new JTextField("?");
            this.jtf_c = new JTextField("?");
            this.jtf_b.setBorder(new EtchedBorder());
            this.jtf_c.setBorder(new EtchedBorder());
            this.jp.add(this.jl_a,"0,0");
            //this.jp.add(new JLabel(" "));
            this.jp.add(this.jtf_b,"2,0");
            this.jp.add(this.jtf_c,"1,0");
            this.add(this.jp,"0,0");

            // add progress bar, if progress is available:
            if(p instanceof AbstractHyperspaceProcess.HasProgress) {
                jpb_a = new JProgressBar();
                jpb_a.setMinimum(0);jpb_a.setMaximum(100);
                jpb_a.setValue(0);
                jpb_a.setStringPainted(true);
                //jpb_a.setForeground(Color.black);
                jpb_a.setBackground(Color.orange);
                //jpb_a.setString("10%");
                this.add(jpb_a,"1,0");
            }

        }

        public void updatePanelGUI() {
            this.jl_a.setText(p.getName());
            this.jtf_b.setText(p.getProcessStatus().name());
            this.jtf_c.setText(p.getProcessStatusMessage());
            if(jpb_a!=null && p instanceof AbstractHyperspaceProcess.HasProgress) {
                int value = (int) (((AbstractHyperspaceProcess.HasProgress)p).getProgress() * 100.0 );
                jpb_a.setValue( value );
                jpb_a.setString(""+value+"%");
                if(value>=100) { jpb_a.setBackground(Color.green); }
                else { jpb_a.setBackground(Color.orange); }
            }


        }

    }

}
