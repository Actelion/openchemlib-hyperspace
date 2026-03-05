package com.idorsia.research.synthonbrowser;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import com.actelion.research.gui.table.ChemistryCellRenderer;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.RawSynthonDownsampler;
import com.idorsia.research.chem.hyperspace.downsampling.SkelSpheresKCentersRawDownsampler;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingRequest;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonSetDownsamplingResult;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSet;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableRowSorter;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
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
    private final DecimalFormat similarityFormat = new DecimalFormat("0.000");

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
        JButton downsampleSetButton = new JButton("Downsample Set...");

        openButton.addActionListener(e -> chooseAndOpen());
        reloadButton.addActionListener(e -> {
            if (currentPath != null) {
                loadSpace(currentPath);
            }
        });
        expandAllButton.addActionListener(e -> expandAll(true));
        collapseAllButton.addActionListener(e -> expandAll(false));
        downsampleSetButton.addActionListener(e -> runDownsamplingForSelectedSet());

        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        panel.add(openButton);
        panel.add(reloadButton);
        panel.add(expandAllButton);
        panel.add(collapseAllButton);
        panel.add(downsampleSetButton);
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

    private SynthonTreeBuilder.NodeData nodeData(DefaultMutableTreeNode node) {
        if (node == null) {
            return null;
        }
        Object userObject = node.getUserObject();
        if (userObject instanceof SynthonTreeBuilder.NodeData nodeData) {
            return nodeData;
        }
        return null;
    }

    private SelectedSet selectedSetFromNode(DefaultMutableTreeNode node) {
        SynthonTreeBuilder.NodeData nodeData = nodeData(node);
        if (nodeData == null || nodeData.type() != SynthonTreeBuilder.NodeType.SET || nodeData.setIndex() == null) {
            return null;
        }
        return new SelectedSet(nodeData.reactionId(), nodeData.setIndex());
    }

    private void runDownsamplingForSelectedSet() {
        runDownsamplingForNode(selectedNode());
    }

    private void runDownsamplingForNode(DefaultMutableTreeNode node) {
        if (currentSpace == null) {
            return;
        }
        SelectedSet selectedSet = selectedSetFromNode(node);
        if (selectedSet == null) {
            JOptionPane.showMessageDialog(this,
                    "Please select a specific synthon set node (Reaction -> Set) first.",
                    "No set selected",
                    JOptionPane.INFORMATION_MESSAGE);
            return;
        }

        List<RawSynthon> sourceSet = getSourceSet(selectedSet);
        if (sourceSet.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "Selected set does not contain synthons.",
                    "Empty set",
                    JOptionPane.INFORMATION_MESSAGE);
            return;
        }

        DownsamplingConfig config = askDownsamplingConfig(selectedSet, sourceSet.size());
        if (config == null) {
            return;
        }

        setBusyState(String.format("Downsampling %s / set %d ...", selectedSet.reactionId(), selectedSet.setIndex()));
        SwingWorker<SynthonSetDownsamplingResult, Void> worker = new SwingWorker<>() {
            @Override
            protected SynthonSetDownsamplingResult doInBackground() {
                SynthonDownsamplingRequest request = SynthonDownsamplingRequest.builder()
                        .withMaxCenters(config.maxCenters())
                        .withSizeCapScale(config.sizeCapScale())
                        .withSizeCapOffset(config.sizeCapOffset())
                        .withMinSimilarity(config.minSimilarity())
                        .withRandomSeed(config.seed())
                        .enforceConnectorEquivalence(config.enforceConnectorEquivalence())
                        .includeClusterMembers(config.includeClusterMembers())
                        .build();
                RawSynthonSet set = new RawSynthonSet(selectedSet.reactionId(), selectedSet.setIndex(), sourceSet);
                RawSynthonDownsampler downsampler = new SkelSpheresKCentersRawDownsampler();
                return downsampler.downsample(set, request);
            }

            @Override
            protected void done() {
                try {
                    SynthonSetDownsamplingResult result = get();
                    statusLabel.setText(String.format(
                            "Downsampling done for %s / set %d: %d -> %d representatives",
                            selectedSet.reactionId(),
                            selectedSet.setIndex(),
                            result.getOriginalSize(),
                            result.getRetainedSize()
                    ));
                    showDownsamplingResultDialog(selectedSet, sourceSet, result);
                } catch (Exception ex) {
                    showError("Could not downsample selected set", ex);
                } finally {
                    setCursor(Cursor.getDefaultCursor());
                }
            }
        };
        worker.execute();
    }

    private List<RawSynthon> getSourceSet(SelectedSet set) {
        RawSynthonSpace.ReactionData reactionData = currentSpace.getReactions().get(set.reactionId());
        if (reactionData == null) {
            return List.of();
        }
        List<RawSynthon> synthons = reactionData.getRawFragmentSets().get(set.setIndex());
        if (synthons == null) {
            return List.of();
        }
        return synthons;
    }

    private DownsamplingConfig askDownsamplingConfig(SelectedSet selectedSet, int sourceSetSize) {
        JTextField maxCentersField = new JTextField("8", 8);
        JTextField sizeCapScaleField = new JTextField("0.0", 8);
        JTextField sizeCapOffsetField = new JTextField("0.0", 8);
        JTextField minSimilarityField = new JTextField("0.75", 8);
        JTextField seedField = new JTextField("13", 8);
        JCheckBox connectorEquivalence = new JCheckBox("Enforce connector equivalence", true);
        JCheckBox includeMembers = new JCheckBox("Capture full cluster members (for expansion view)", true);

        JPanel panel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 6, 4, 6);
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;

        panel.add(new JLabel(String.format("Reaction: %s   Set: %d   Size: %d",
                selectedSet.reactionId(),
                selectedSet.setIndex(),
                sourceSetSize)), gbc);

        gbc.gridy++;
        panel.add(new JLabel("maxCenters (0 = unlimited)"), gbc);
        gbc.gridx = 1;
        panel.add(maxCentersField, gbc);

        gbc.gridx = 0;
        gbc.gridy++;
        panel.add(new JLabel("sizeCapScale"), gbc);
        gbc.gridx = 1;
        panel.add(sizeCapScaleField, gbc);

        gbc.gridx = 0;
        gbc.gridy++;
        panel.add(new JLabel("sizeCapOffset"), gbc);
        gbc.gridx = 1;
        panel.add(sizeCapOffsetField, gbc);

        gbc.gridx = 0;
        gbc.gridy++;
        panel.add(new JLabel("minSimilarity [0..1]"), gbc);
        gbc.gridx = 1;
        panel.add(minSimilarityField, gbc);

        gbc.gridx = 0;
        gbc.gridy++;
        panel.add(new JLabel("seed"), gbc);
        gbc.gridx = 1;
        panel.add(seedField, gbc);

        gbc.gridx = 0;
        gbc.gridy++;
        gbc.gridwidth = 2;
        panel.add(connectorEquivalence, gbc);
        gbc.gridy++;
        panel.add(includeMembers, gbc);

        while (true) {
            int rc = JOptionPane.showConfirmDialog(
                    this,
                    panel,
                    "Downsample Set - SkelSpheresKCentersRaw",
                    JOptionPane.OK_CANCEL_OPTION,
                    JOptionPane.PLAIN_MESSAGE
            );
            if (rc != JOptionPane.OK_OPTION) {
                return null;
            }

            try {
                int maxCenters = Integer.parseInt(maxCentersField.getText().trim());
                double sizeCapScale = Double.parseDouble(sizeCapScaleField.getText().trim());
                double sizeCapOffset = Double.parseDouble(sizeCapOffsetField.getText().trim());
                double minSimilarity = Double.parseDouble(minSimilarityField.getText().trim());
                long seed = Long.parseLong(seedField.getText().trim());

                if (maxCenters < 0) {
                    throw new IllegalArgumentException("maxCenters must be >= 0");
                }
                if (sizeCapScale < 0.0) {
                    throw new IllegalArgumentException("sizeCapScale must be >= 0");
                }
                if (minSimilarity < 0.0 || minSimilarity > 1.0) {
                    throw new IllegalArgumentException("minSimilarity must be in [0, 1]");
                }

                return new DownsamplingConfig(
                        maxCenters,
                        sizeCapScale,
                        sizeCapOffset,
                        minSimilarity,
                        seed,
                        connectorEquivalence.isSelected(),
                        includeMembers.isSelected()
                );
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(
                        this,
                        "Invalid settings: " + ex.getMessage(),
                        "Validation Error",
                        JOptionPane.ERROR_MESSAGE
                );
            }
        }
    }

    private void showDownsamplingResultDialog(SelectedSet selectedSet,
                                              List<RawSynthon> sourceSet,
                                              SynthonSetDownsamplingResult result) {
        Map<String, RawSynthon> byFragmentId = new HashMap<>();
        for (RawSynthon synthon : sourceSet) {
            byFragmentId.put(synthon.getFragmentId(), synthon);
        }

        List<RepresentativeRow> representativeRows = new ArrayList<>();
        for (SynthonSetDownsamplingResult.ClusterInfo cluster : result.getClusterInfos()) {
            RawSynthon representative = resolveRawSynthon(cluster.getRepresentative(), byFragmentId);
            representativeRows.add(new RepresentativeRow(
                    representative.getIdcode(),
                    representative.getFragmentId(),
                    SynthonTableModel.formatConnectors(representative.getConnectors()),
                    cluster.getMembers(),
                    cluster.getMinSimilarityWithinCluster(),
                    cluster.getMemberAssignments().size(),
                    cluster
            ));
        }
        representativeRows.sort(Comparator.comparing(RepresentativeRow::fragmentId, Comparator.nullsFirst(Comparator.naturalOrder())));

        RepresentativeTableModel representativesModel = new RepresentativeTableModel(representativeRows);
        ClusterMemberTableModel memberTableModel = new ClusterMemberTableModel();

        JTable representativesTable = new JTable(representativesModel);
        representativesTable.setRowHeight(90);
        representativesTable.getColumnModel().getColumn(0).setCellRenderer(new ChemistryCellRenderer());
        representativesTable.getColumnModel().getColumn(0).setPreferredWidth(220);
        representativesTable.getColumnModel().getColumn(2).setPreferredWidth(120);

        JTable membersTable = new JTable(memberTableModel);
        membersTable.setRowHeight(90);
        membersTable.getColumnModel().getColumn(0).setCellRenderer(new ChemistryCellRenderer());
        membersTable.getColumnModel().getColumn(0).setPreferredWidth(220);
        membersTable.getColumnModel().getColumn(2).setPreferredWidth(120);
        membersTable.getColumnModel().getColumn(3).setPreferredWidth(110);

        JLabel membersHeader = new JLabel("Cluster Members");
        representativesTable.getSelectionModel().addListSelectionListener(e -> {
            if (e.getValueIsAdjusting()) {
                return;
            }
            int row = representativesTable.getSelectedRow();
            if (row < 0) {
                memberTableModel.setRows(List.of());
                membersHeader.setText("Cluster Members");
                return;
            }
            RepresentativeRow selected = representativesModel.getRow(row);
            List<MemberRow> memberRows = toMemberRows(selected.clusterInfo(), byFragmentId);
            memberTableModel.setRows(memberRows);
            if (selected.clusterInfo().hasMemberAssignments()) {
                membersHeader.setText(String.format("Cluster Members (%d captured)", memberRows.size()));
            } else {
                membersHeader.setText("Cluster Members (not captured for this run)");
            }
        });

        JLabel summaryLabel = new JLabel(String.format(
                "%s / Set %d: %d -> %d representatives (discarded=%d)",
                selectedSet.reactionId(),
                selectedSet.setIndex(),
                result.getOriginalSize(),
                result.getRetainedSize(),
                result.getDiscardedSize()
        ));

        JPanel top = new JPanel(new BorderLayout());
        top.add(summaryLabel, BorderLayout.NORTH);
        top.add(new JScrollPane(representativesTable), BorderLayout.CENTER);

        JPanel bottom = new JPanel(new BorderLayout());
        bottom.add(membersHeader, BorderLayout.NORTH);
        bottom.add(new JScrollPane(membersTable), BorderLayout.CENTER);

        JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT, top, bottom);
        split.setResizeWeight(0.55);

        JDialog dialog = new JDialog(this, "Downsampling Result", true);
        dialog.setLayout(new BorderLayout());
        dialog.add(split, BorderLayout.CENTER);
        dialog.setSize(1150, 760);
        dialog.setLocationRelativeTo(this);
        if (!representativeRows.isEmpty()) {
            representativesTable.setRowSelectionInterval(0, 0);
        }
        dialog.setVisible(true);
    }

    private List<MemberRow> toMemberRows(SynthonSetDownsamplingResult.ClusterInfo clusterInfo,
                                         Map<String, RawSynthon> byFragmentId) {
        if (!clusterInfo.hasMemberAssignments()) {
            return List.of();
        }
        List<MemberRow> out = new ArrayList<>();
        for (SynthonSetDownsamplingResult.ClusterMember member : clusterInfo.getMemberAssignments()) {
            RawSynthon raw = resolveRawSynthon(member.getFragId(), byFragmentId);
            out.add(new MemberRow(
                    raw.getIdcode(),
                    raw.getFragmentId(),
                    SynthonTableModel.formatConnectors(raw.getConnectors()),
                    member.getSimilarityToRepresentative()
            ));
        }
        out.sort(Comparator.comparingDouble(MemberRow::similarityToRepresentative).reversed());
        return out;
    }

    private RawSynthon resolveRawSynthon(SynthonSpace.FragId fragId, Map<String, RawSynthon> byFragmentId) {
        RawSynthon fromMap = byFragmentId.get(fragId.fragment_id);
        if (fromMap != null) {
            return fromMap;
        }
        return RawSynthon.fromFragId(fragId);
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
        SynthonTreeBuilder.NodeData nodeData = nodeData(node);
        if (nodeData != null && nodeData.type() == SynthonTreeBuilder.NodeType.SET) {
            JMenuItem downsampleSet = new JMenuItem("Downsample set...");
            downsampleSet.addActionListener(a -> runDownsamplingForNode(node));
            menu.add(downsampleSet);
            menu.addSeparator();
        }
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

    private record SelectedSet(String reactionId, int setIndex) {}

    private record DownsamplingConfig(int maxCenters,
                                      double sizeCapScale,
                                      double sizeCapOffset,
                                      double minSimilarity,
                                      long seed,
                                      boolean enforceConnectorEquivalence,
                                      boolean includeClusterMembers) {}

    private record RepresentativeRow(String idcode,
                                     String fragmentId,
                                     String connectors,
                                     int members,
                                     double minSimilarityWithinCluster,
                                     int capturedAssignments,
                                     SynthonSetDownsamplingResult.ClusterInfo clusterInfo) {}

    private record MemberRow(String idcode,
                             String fragmentId,
                             String connectors,
                             double similarityToRepresentative) {}

    private final class RepresentativeTableModel extends AbstractTableModel {

        private final String[] columns = {
                "Representative[idcode]",
                "Fragment ID",
                "Connectors",
                "Cluster Members",
                "Min Similarity",
                "Captured Assignments"
        };
        private final List<RepresentativeRow> rows;

        private RepresentativeTableModel(List<RepresentativeRow> rows) {
            this.rows = new ArrayList<>(rows);
        }

        @Override
        public int getRowCount() {
            return rows.size();
        }

        @Override
        public int getColumnCount() {
            return columns.length;
        }

        @Override
        public String getColumnName(int column) {
            return columns[column];
        }

        @Override
        public Class<?> getColumnClass(int columnIndex) {
            return switch (columnIndex) {
                case 3, 5 -> Integer.class;
                case 4 -> String.class;
                default -> String.class;
            };
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            RepresentativeRow row = rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.idcode();
                case 1 -> row.fragmentId();
                case 2 -> row.connectors();
                case 3 -> row.members();
                case 4 -> similarityFormat.format(row.minSimilarityWithinCluster());
                case 5 -> row.capturedAssignments();
                default -> "";
            };
        }

        public RepresentativeRow getRow(int rowIndex) {
            return rows.get(rowIndex);
        }
    }

    private final class ClusterMemberTableModel extends AbstractTableModel {

        private final String[] columns = {
                "Member[idcode]",
                "Fragment ID",
                "Connectors",
                "Similarity to Rep"
        };
        private final List<MemberRow> rows = new ArrayList<>();

        @Override
        public int getRowCount() {
            return rows.size();
        }

        @Override
        public int getColumnCount() {
            return columns.length;
        }

        @Override
        public String getColumnName(int column) {
            return columns[column];
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            MemberRow row = rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.idcode();
                case 1 -> row.fragmentId();
                case 2 -> row.connectors();
                case 3 -> similarityFormat.format(row.similarityToRepresentative());
                default -> "";
            };
        }

        public void setRows(List<MemberRow> memberRows) {
            rows.clear();
            rows.addAll(memberRows);
            fireTableDataChanged();
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
