package com.idorsia.research.synthonbrowser;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import com.actelion.research.gui.table.ChemistryCellRenderer;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.TableRowSorter;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

public class SynthonSpaceBrowserFrame extends JFrame {

    private final JTree tree = new JTree(new DefaultMutableTreeNode("No synthon space loaded"));
    private final SynthonTableModel tableModel = new SynthonTableModel();
    private final JTable synthonTable = new JTable(tableModel);
    private final TableRowSorter<SynthonTableModel> tableSorter = new TableRowSorter<>(tableModel);
    private final JTextField filterField = new JTextField(28);
    private final JCheckBox downsampledSetsCheckBox = new JCheckBox("Use downsampled sets");
    private final JLabel statusLabel = new JLabel("Load a RawSynthonSpace to start.");
    private final JPanel previewPanel = new JPanel(new BorderLayout());

    private RawSynthonSpace currentSpace;
    private Path currentPath;

    public SynthonSpaceBrowserFrame() {
        super("Synthon Space Browser");
        initUI();
    }

    public void openPath(Path path) {
        loadSpace(path);
    }

    private void initUI() {
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setLayout(new BorderLayout());

        add(buildToolbar(), BorderLayout.NORTH);
        add(buildMainContent(), BorderLayout.CENTER);
        add(statusLabel, BorderLayout.SOUTH);

        tree.addTreeSelectionListener(e -> refreshTableFromTreeSelection());
        tree.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                maybeShowTreeContextMenu(e);
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                maybeShowTreeContextMenu(e);
            }
        });
        synthonTable.getSelectionModel().addListSelectionListener(e -> {
            if (!e.getValueIsAdjusting()) {
                updatePreviewFromSelection();
            }
        });
        filterField.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                applyFilter();
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                applyFilter();
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                applyFilter();
            }
        });
        downsampledSetsCheckBox.addActionListener(e -> rebuildTreeModel());

        setPreferredSize(new Dimension(1300, 820));
        pack();
        setLocationRelativeTo(null);
    }

    private JComponent buildToolbar() {
        JButton openButton = new JButton("Open...");
        JButton reloadButton = new JButton("Reload");
        JButton expandAllButton = new JButton("Expand All");
        JButton collapseAllButton = new JButton("Collapse All");

        openButton.addActionListener(e -> chooseAndOpen());
        reloadButton.addActionListener(e -> {
            if (currentPath != null) {
                loadSpace(currentPath);
            }
        });
        expandAllButton.addActionListener(e -> expandAll(true));
        collapseAllButton.addActionListener(e -> expandAll(false));

        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        panel.add(openButton);
        panel.add(reloadButton);
        panel.add(expandAllButton);
        panel.add(collapseAllButton);
        panel.add(downsampledSetsCheckBox);
        panel.add(new JLabel("Filter:"));
        panel.add(filterField);
        return panel;
    }

    private JComponent buildMainContent() {
        tree.setRootVisible(true);
        JScrollPane treeScroll = new JScrollPane(tree);
        treeScroll.setPreferredSize(new Dimension(380, 600));

        synthonTable.setRowHeight(96);
        synthonTable.setAutoCreateRowSorter(false);
        synthonTable.setRowSorter(tableSorter);
        synthonTable.getColumnModel().getColumn(0).setCellRenderer(new ChemistryCellRenderer());
        synthonTable.getColumnModel().getColumn(0).setPreferredWidth(220);
        synthonTable.getColumnModel().getColumn(3).setPreferredWidth(300);

        JScrollPane tableScroll = new JScrollPane(synthonTable);

        previewPanel.setBorder(BorderFactory.createTitledBorder("Preview"));
        previewPanel.setPreferredSize(new Dimension(400, 200));
        setPreviewMolecule(new StereoMolecule());

        JSplitPane rightSplit = new JSplitPane(JSplitPane.VERTICAL_SPLIT, tableScroll, previewPanel);
        rightSplit.setResizeWeight(0.8);

        JSplitPane mainSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treeScroll, rightSplit);
        mainSplit.setResizeWeight(0.27);
        return mainSplit;
    }

    private void chooseAndOpen() {
        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle("Open RawSynthonSpace file");
        int rc = chooser.showOpenDialog(this);
        if (rc == JFileChooser.APPROVE_OPTION) {
            loadSpace(chooser.getSelectedFile().toPath());
        }
    }

    private void loadSpace(Path path) {
        Objects.requireNonNull(path, "path");
        setBusyState("Loading: " + path.toAbsolutePath());
        SwingWorker<RawSynthonSpace, Void> worker = new SwingWorker<>() {
            @Override
            protected RawSynthonSpace doInBackground() throws Exception {
                return RawSynthonSpaceIO.read(path);
            }

            @Override
            protected void done() {
                try {
                    currentSpace = get();
                    currentPath = path;
                    boolean hasDownsampled = currentSpace.hasDownsampledSets();
                    downsampledSetsCheckBox.setEnabled(hasDownsampled);
                    if (!hasDownsampled) {
                        downsampledSetsCheckBox.setSelected(false);
                    }
                    rebuildTreeModel();
                    statusLabel.setText(String.format(
                            "Loaded %s (reactions=%d, mode=%s)",
                            path.toAbsolutePath(),
                            currentSpace.getReactions().size(),
                            downsampledSetsCheckBox.isSelected() ? "downsampled" : "full"
                    ));
                } catch (Exception ex) {
                    showError("Could not load synthon space", ex);
                } finally {
                    setCursor(Cursor.getDefaultCursor());
                }
            }
        };
        worker.execute();
    }

    private void rebuildTreeModel() {
        if (currentSpace == null) {
            return;
        }
        boolean useDownsampled = downsampledSetsCheckBox.isSelected();
        DefaultMutableTreeNode root = SynthonTreeBuilder.buildTree(currentSpace, useDownsampled);
        tree.setModel(new DefaultTreeModel(root));
        expandFirstTwoLevels();
        tree.setSelectionPath(new TreePath(root.getPath()));
    }

    private void refreshTableFromTreeSelection() {
        if (currentSpace == null) {
            tableModel.setRows(List.of());
            return;
        }
        DefaultMutableTreeNode node = selectedNode();
        List<SynthonTableModel.Row> rows = collectRowsForNode(node);
        sortRows(rows);
        tableModel.setRows(rows);
        applyFilter();
        if (!rows.isEmpty()) {
            synthonTable.setRowSelectionInterval(0, 0);
        } else {
            setPreviewMolecule(new StereoMolecule());
        }
        statusLabel.setText(String.format(
                "Rows: %d (view: %d)",
                rows.size(),
                synthonTable.getRowCount()
        ));
    }

    private void sortRows(List<SynthonTableModel.Row> rows) {
        rows.sort(
                Comparator.comparing(SynthonTableModel.Row::reactionId)
                        .thenComparingInt(SynthonTableModel.Row::setIndex)
                        .thenComparing(row -> row.synthon().getFragmentId(), Comparator.nullsFirst(Comparator.naturalOrder()))
        );
    }

    private List<SynthonTableModel.Row> collectRowsForNode(DefaultMutableTreeNode node) {
        if (node == null) {
            return collectAllRows();
        }
        Object userObject = node.getUserObject();
        if (!(userObject instanceof SynthonTreeBuilder.NodeData nodeData)) {
            return collectAllRows();
        }
        return switch (nodeData.type()) {
            case ROOT -> collectAllRows();
            case REACTION -> collectRowsForReaction(nodeData.reactionId());
            case SET -> collectRowsForSet(nodeData.reactionId(), nodeData.setIndex());
        };
    }

    private List<SynthonTableModel.Row> collectAllRows() {
        List<SynthonTableModel.Row> out = new ArrayList<>();
        for (Map.Entry<String, RawSynthonSpace.ReactionData> reactionEntry : currentSpace.getReactions().entrySet()) {
            out.addAll(rowsForReactionSets(reactionEntry.getKey(), reactionEntry.getValue()));
        }
        return out;
    }

    private List<SynthonTableModel.Row> collectRowsForReaction(String reactionId) {
        RawSynthonSpace.ReactionData rd = currentSpace.getReactions().get(reactionId);
        if (rd == null) {
            return List.of();
        }
        return rowsForReactionSets(reactionId, rd);
    }

    private List<SynthonTableModel.Row> collectRowsForSet(String reactionId, Integer setIndex) {
        if (setIndex == null) {
            return List.of();
        }
        RawSynthonSpace.ReactionData rd = currentSpace.getReactions().get(reactionId);
        if (rd == null) {
            return List.of();
        }
        Map<Integer, List<RawSynthon>> sets = SynthonTreeBuilder.activeSets(rd, downsampledSetsCheckBox.isSelected());
        List<RawSynthon> synthons = sets.get(setIndex);
        if (synthons == null || synthons.isEmpty()) {
            return List.of();
        }
        List<SynthonTableModel.Row> out = new ArrayList<>(synthons.size());
        for (RawSynthon synthon : synthons) {
            out.add(new SynthonTableModel.Row(reactionId, setIndex, synthon));
        }
        return out;
    }

    private List<SynthonTableModel.Row> rowsForReactionSets(String reactionId, RawSynthonSpace.ReactionData rd) {
        Map<Integer, List<RawSynthon>> sets = SynthonTreeBuilder.activeSets(rd, downsampledSetsCheckBox.isSelected());
        List<Integer> setIds = new ArrayList<>(sets.keySet());
        setIds.sort(Comparator.naturalOrder());
        List<SynthonTableModel.Row> out = new ArrayList<>();
        for (Integer setId : setIds) {
            List<RawSynthon> synthons = sets.get(setId);
            if (synthons == null) {
                continue;
            }
            for (RawSynthon synthon : synthons) {
                out.add(new SynthonTableModel.Row(reactionId, setId, synthon));
            }
        }
        return out;
    }

    private void maybeShowTreeContextMenu(MouseEvent e) {
        if (!e.isPopupTrigger()) {
            return;
        }
        if (currentSpace == null) {
            return;
        }
        TreePath path = tree.getPathForLocation(e.getX(), e.getY());
        if (path == null) {
            return;
        }
        tree.setSelectionPath(path);
        Object nodeObj = path.getLastPathComponent();
        if (!(nodeObj instanceof DefaultMutableTreeNode node)) {
            return;
        }

        JPopupMenu menu = new JPopupMenu();
        JMenuItem exportTsv = new JMenuItem("Export to .tsv");
        exportTsv.addActionListener(a -> exportNodeToTsv(node));
        menu.add(exportTsv);
        menu.show(tree, e.getX(), e.getY());
    }

    private void exportNodeToTsv(DefaultMutableTreeNode node) {
        List<SynthonTableModel.Row> rows = collectRowsForNode(node);
        sortRows(rows);

        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle("Export synthons to TSV");
        chooser.setFileFilter(new FileNameExtensionFilter("TSV files", "tsv"));
        chooser.setSelectedFile(suggestExportFile(node));
        int rc = chooser.showSaveDialog(this);
        if (rc != JFileChooser.APPROVE_OPTION) {
            return;
        }

        Path out = chooser.getSelectedFile().toPath();
        if (!out.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".tsv")) {
            out = out.resolveSibling(out.getFileName().toString() + ".tsv");
        }

        Path finalOut = out;
        setBusyState("Exporting TSV: " + finalOut.toAbsolutePath());
        SwingWorker<Void, Void> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() throws Exception {
                writeRowsToTsv(rows, finalOut);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                    statusLabel.setText(String.format(
                            "Exported %d rows to %s",
                            rows.size(),
                            finalOut.toAbsolutePath()
                    ));
                } catch (Exception ex) {
                    showError("Could not export TSV", ex);
                } finally {
                    setCursor(Cursor.getDefaultCursor());
                }
            }
        };
        worker.execute();
    }

    private java.io.File suggestExportFile(DefaultMutableTreeNode node) {
        String base = "synthons";
        Object userObject = node.getUserObject();
        if (userObject instanceof SynthonTreeBuilder.NodeData nodeData) {
            base = switch (nodeData.type()) {
                case ROOT -> "synthons_all";
                case REACTION -> "synthons_" + safeFileToken(nodeData.reactionId());
                case SET -> "synthons_" + safeFileToken(nodeData.reactionId()) + "_set_" + nodeData.setIndex();
            };
        }
        return new java.io.File(base + ".tsv");
    }

    private String safeFileToken(String value) {
        if (value == null || value.isBlank()) {
            return "na";
        }
        return value.replaceAll("[^a-zA-Z0-9._-]+", "_");
    }

    private void writeRowsToTsv(List<SynthonTableModel.Row> rows, Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            for (int c = 0; c < tableModel.getColumnCount(); c++) {
                if (c > 0) {
                    writer.write('\t');
                }
                writer.write(escapeTsv(tableModel.getColumnName(c)));
            }
            writer.newLine();

            for (SynthonTableModel.Row row : rows) {
                for (int c = 0; c < tableModel.getColumnCount(); c++) {
                    if (c > 0) {
                        writer.write('\t');
                    }
                    Object value = tableModel.getValueAt(row, c);
                    writer.write(escapeTsv(value == null ? "" : String.valueOf(value)));
                }
                writer.newLine();
            }
        }
    }

    private String escapeTsv(String value) {
        if (value == null) {
            return "";
        }
        return value
                .replace('\t', ' ')
                .replace('\r', ' ')
                .replace('\n', ' ');
    }

    private void applyFilter() {
        String text = filterField.getText();
        if (text == null || text.isBlank()) {
            tableSorter.setRowFilter(null);
            return;
        }
        String needle = text.toLowerCase(Locale.ROOT);
        tableSorter.setRowFilter(new RowFilter<>() {
            @Override
            public boolean include(Entry<? extends SynthonTableModel, ? extends Integer> entry) {
                for (int c = 1; c < tableModel.getColumnCount(); c++) {
                    Object v = entry.getValue(c);
                    if (v != null && v.toString().toLowerCase(Locale.ROOT).contains(needle)) {
                        return true;
                    }
                }
                return false;
            }
        });
    }

    private void updatePreviewFromSelection() {
        int viewRow = synthonTable.getSelectedRow();
        if (viewRow < 0) {
            setPreviewMolecule(new StereoMolecule());
            return;
        }
        int modelRow = synthonTable.convertRowIndexToModel(viewRow);
        SynthonTableModel.Row row = tableModel.getRow(modelRow);
        StereoMolecule molecule = parseIdcode(row.synthon().getIdcode());
        setPreviewMolecule(molecule);
    }

    private StereoMolecule parseIdcode(String idcode) {
        StereoMolecule molecule = new StereoMolecule();
        if (idcode == null || idcode.isBlank()) {
            return molecule;
        }
        try {
            IDCodeParser parser = new IDCodeParser();
            parser.parse(molecule, idcode);
            molecule.ensureHelperArrays(StereoMolecule.cHelperCIP);
        } catch (Exception ignored) {
            return new StereoMolecule();
        }
        return molecule;
    }

    private void setPreviewMolecule(StereoMolecule molecule) {
        previewPanel.removeAll();
        previewPanel.add(new JStructureView(molecule), BorderLayout.CENTER);
        previewPanel.revalidate();
        previewPanel.repaint();
    }

    private DefaultMutableTreeNode selectedNode() {
        TreePath path = tree.getSelectionPath();
        if (path == null) {
            return null;
        }
        Object node = path.getLastPathComponent();
        return node instanceof DefaultMutableTreeNode dmtn ? dmtn : null;
    }

    private void expandFirstTwoLevels() {
        if (tree.getRowCount() == 0) {
            return;
        }
        tree.expandRow(0);
        for (int i = 1; i < tree.getRowCount(); i++) {
            TreePath path = tree.getPathForRow(i);
            if (path != null && path.getPathCount() <= 3) {
                tree.expandRow(i);
            }
        }
    }

    private void expandAll(boolean expand) {
        int row = 0;
        while (row < tree.getRowCount()) {
            if (expand) {
                tree.expandRow(row);
            } else if (row > 0) {
                tree.collapseRow(row);
            }
            row++;
        }
    }

    private void setBusyState(String message) {
        statusLabel.setText(message);
        setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    }

    private void showError(String prefix, Exception ex) {
        String msg = ex.getMessage() == null ? ex.toString() : ex.getMessage();
        statusLabel.setText(prefix + ": " + msg);
        JOptionPane.showMessageDialog(this, prefix + "\n" + msg, "Error", JOptionPane.ERROR_MESSAGE);
    }
}
