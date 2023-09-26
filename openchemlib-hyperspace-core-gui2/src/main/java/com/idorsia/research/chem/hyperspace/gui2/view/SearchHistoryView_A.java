package com.idorsia.research.chem.hyperspace.gui2.view;

import com.idorsia.research.chem.hyperspace.gui2.model.EnamineHyperspaceModel_A;
import com.idorsia.research.chem.hyperspace.gui2.model.ProcessTableModel;
import com.idorsia.research.chem.hyperspace.gui2.model.RealTimeExpandingSearchResultModel;
import com.idorsia.research.chem.hyperspace.gui2.task.LongRunningTask;
import com.idorsia.research.chem.hyperspace.gui2.task.SubstructureSearchTask;

import javax.swing.*;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

/**
 * Is a process table that just shows the past searches
 */
public class SearchHistoryView_A extends JPanel {

        private EnamineHyperspaceModel_A hyperspaceModel;
        private ProcessTableModel model;

        private JPanel panelTop;
        private JTable table;

    public SearchHistoryView_A(EnamineHyperspaceModel_A hyperspaceModel, ProcessTableModel model) {
        this.hyperspaceModel = hyperspaceModel;
        this.model = model;
        reinit();
    }

        private void reinit() {
            this.removeAll();
            setLayout(new BorderLayout());

            // Create the table with three columns: Name, Status, Progress
            String[] columnNames = {"Query", "Status", "Progress"};
            table = new JTable(model.getPastSearchesTableModel());
            JScrollPane scrollPane = new JScrollPane(table);
            add(scrollPane, BorderLayout.CENTER);

            panelTop = new JPanel(); panelTop.setLayout(new FlowLayout(FlowLayout.LEFT));
            panelTop.add(new JLabel("Search History"));
            this.add(panelTop,BorderLayout.NORTH);

            // configure table:
            table.setRowHeight(80);
            FragmentHighlightChemistryCellRenderer ccr = new FragmentHighlightChemistryCellRenderer(new Dimension(80,120),null);
            this.table.getColumnModel().getColumn(0).setCellRenderer(ccr);

            model.getPastSearchesTableModel().addTableModelListener(new TableModelListener() {
                @Override
                public void tableChanged(TableModelEvent e) {
                    SwingUtilities.updateComponentTreeUI(table);
                }
            });


            // init mouse context menu for jtable
            initMouseContextMenu();

        }

        private void initMouseContextMenu() {
            JPopupMenu popupMenu = new JPopupMenu();
            JMenuItem showSearchResults = new JMenuItem("Show results");
            JMenuItem stopSearch        = new JMenuItem("Stop search");
            showSearchResults.addActionListener( (e) -> {
                int selectedRow = table.getSelectedRow();
                if (selectedRow != -1) {
                    hyperspaceModel.setPresentedSearchResult( new RealTimeExpandingSearchResultModel( model.getPastSearches().get(selectedRow).getResultsModel() , 32000 ) );
                }
            } );
            stopSearch.addActionListener( (e) -> {
                int selectedRow = table.getSelectedRow();
                if (selectedRow != -1) {
                    SubstructureSearchTask ti = model.getPastSearches().get(selectedRow);
                    ti.cancel(true);
                }
            } );
            popupMenu.add(showSearchResults);
            popupMenu.add(stopSearch);

            table.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseReleased(MouseEvent e) {
                    if (SwingUtilities.isRightMouseButton(e) && table.rowAtPoint(e.getPoint()) != -1) {
                        table.setRowSelectionInterval(table.rowAtPoint(e.getPoint()), table.rowAtPoint(e.getPoint()));
                        popupMenu.show(e.getComponent(), e.getX(), e.getY());
                    }
                }
            });
        }


        public static void main(String[] args) {
            SwingUtilities.invokeLater(() -> {
                JFrame frame = new JFrame("Process Table Example");
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.setSize(400, 300);

                ProcessTableModel pmodel = new ProcessTableModel();
                ProcessTableView panel = new ProcessTableView(pmodel);
                frame.add(panel);

                // Create and add some tasks to the panel
                LongRunningTask task1 = new LongRunningTask("Task 1");
                LongRunningTask task2 = new LongRunningTask("Task 2");
                pmodel.addTask(task1);
                pmodel.addTask(task2);

                task1.execute();

                frame.setVisible(true);
            });
        }
}
