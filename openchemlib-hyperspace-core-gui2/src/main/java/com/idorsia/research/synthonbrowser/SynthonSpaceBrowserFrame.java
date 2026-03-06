package com.idorsia.research.synthonbrowser;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import com.actelion.research.gui.table.ChemistryCellRenderer;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
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
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;

public class SynthonSpaceBrowserFrame extends JFrame {

    private static final String CARD_SYNTHONS = "synthons";
    private static final String CARD_REACTION = "reaction";
    private static final String CARD_SPACE = "space";
    private static final int RANDOM_EXAMPLE_COUNT = 10;
    private static final int RANDOM_SAMPLE_MAX_ATTEMPTS = 300;
    private static final int SPACE_LARGEST_REACTIONS_LIMIT = 20;
    private static final String[] SPACE_BUCKET_LABELS = {
            "<1M",
            "1M-10M",
            "10M-100M",
            "100M-1B",
            "1B-10B",
            "10B-100B",
            ">100B"
    };

    private final JTree tree = new JTree(new DefaultMutableTreeNode("No synthon space loaded"));
    private final SynthonTableModel tableModel = new SynthonTableModel();
    private final JTable synthonTable = new JTable(tableModel);
    private final TableRowSorter<SynthonTableModel> tableSorter = new TableRowSorter<>(tableModel);
    private final JTextField filterField = new JTextField(28);
    private final JCheckBox downsampledSetsCheckBox = new JCheckBox("Use downsampled sets");
    private final JLabel statusLabel = new JLabel("Load a RawSynthonSpace to start.");
    private final JPanel previewPanel = new JPanel(new BorderLayout());
    private final CardLayout rightTopCards = new CardLayout();
    private final JPanel rightTopPanel = new JPanel(rightTopCards);
    private final JLabel reactionSummaryLabel = new JLabel("Select a reaction node to inspect statistics.");
    private final ReactionHistogramTableModel reactionHistogramModel = new ReactionHistogramTableModel();
    private final JTable reactionHistogramTable = new JTable(reactionHistogramModel);
    private final RandomExamplesTableModel randomExamplesTableModel = new RandomExamplesTableModel();
    private final JTable randomExamplesTable = new JTable(randomExamplesTableModel);
    private final JLabel spaceSummaryLabel = new JLabel("Select the space node to inspect global statistics.");
    private final SpaceHistogramTableModel spaceHistogramModel = new SpaceHistogramTableModel();
    private final JTable spaceHistogramTable = new JTable(spaceHistogramModel);
    private final LargestReactionsTableModel largestReactionsTableModel = new LargestReactionsTableModel();
    private final JTable largestReactionsTable = new JTable(largestReactionsTableModel);
    private final Map<String, Integer> nonHydrogenAtomCache = new HashMap<>();
    private final Map<String, ReactionStats> reactionStatsCache = new HashMap<>();
    private final Map<String, SpaceStats> spaceStatsCache = new HashMap<>();
    private final DecimalFormat similarityFormat = new DecimalFormat("0.000");
    private final DecimalFormat percentFormat = new DecimalFormat("0.00");

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
        randomExamplesTable.getSelectionModel().addListSelectionListener(e -> {
            if (!e.getValueIsAdjusting()) {
                updatePreviewFromRandomExampleSelection();
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
        downsampledSetsCheckBox.addActionListener(e -> {
            reactionStatsCache.clear();
            spaceStatsCache.clear();
            rebuildTreeModel();
        });
        rightTopCards.show(rightTopPanel, CARD_SYNTHONS);

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
        rightTopPanel.add(tableScroll, CARD_SYNTHONS);
        rightTopPanel.add(buildReactionDetailsPanel(), CARD_REACTION);
        rightTopPanel.add(buildSpaceDetailsPanel(), CARD_SPACE);

        previewPanel.setBorder(BorderFactory.createTitledBorder("Preview"));
        previewPanel.setPreferredSize(new Dimension(400, 200));
        setPreviewMolecule(new StereoMolecule());

        JSplitPane rightSplit = new JSplitPane(JSplitPane.VERTICAL_SPLIT, rightTopPanel, previewPanel);
        rightSplit.setResizeWeight(0.8);

        JSplitPane mainSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treeScroll, rightSplit);
        mainSplit.setResizeWeight(0.27);
        return mainSplit;
    }

    private JComponent buildReactionDetailsPanel() {
        reactionSummaryLabel.setBorder(BorderFactory.createEmptyBorder(6, 8, 6, 8));

        reactionHistogramTable.setFillsViewportHeight(true);
        reactionHistogramTable.getColumnModel().getColumn(0).setPreferredWidth(120);
        reactionHistogramTable.getColumnModel().getColumn(1).setPreferredWidth(220);
        reactionHistogramTable.getColumnModel().getColumn(2).setPreferredWidth(140);
        reactionHistogramTable.getColumnModel().getColumn(3).setPreferredWidth(120);

        randomExamplesTable.setRowHeight(90);
        randomExamplesTable.setFillsViewportHeight(true);
        randomExamplesTable.getColumnModel().getColumn(0).setCellRenderer(new ChemistryCellRenderer());
        randomExamplesTable.getColumnModel().getColumn(0).setPreferredWidth(240);
        randomExamplesTable.getColumnModel().getColumn(1).setPreferredWidth(120);
        randomExamplesTable.getColumnModel().getColumn(2).setPreferredWidth(260);

        JTabbedPane tabs = new JTabbedPane();
        tabs.addTab("Size Histogram", new JScrollPane(reactionHistogramTable));
        tabs.addTab("Random Compounds", new JScrollPane(randomExamplesTable));

        JPanel panel = new JPanel(new BorderLayout());
        panel.add(reactionSummaryLabel, BorderLayout.NORTH);
        panel.add(tabs, BorderLayout.CENTER);
        return panel;
    }

    private JComponent buildSpaceDetailsPanel() {
        spaceSummaryLabel.setBorder(BorderFactory.createEmptyBorder(6, 8, 6, 8));

        spaceHistogramTable.setFillsViewportHeight(true);
        spaceHistogramTable.getColumnModel().getColumn(0).setPreferredWidth(160);
        spaceHistogramTable.getColumnModel().getColumn(1).setPreferredWidth(140);
        spaceHistogramTable.getColumnModel().getColumn(2).setPreferredWidth(110);

        largestReactionsTable.setFillsViewportHeight(true);
        largestReactionsTable.getColumnModel().getColumn(0).setPreferredWidth(280);
        largestReactionsTable.getColumnModel().getColumn(1).setPreferredWidth(70);
        largestReactionsTable.getColumnModel().getColumn(2).setPreferredWidth(100);
        largestReactionsTable.getColumnModel().getColumn(3).setPreferredWidth(130);

        JTabbedPane tabs = new JTabbedPane();
        tabs.addTab("Reaction Size Buckets", new JScrollPane(spaceHistogramTable));
        tabs.addTab("Largest Reactions", new JScrollPane(largestReactionsTable));

        JPanel panel = new JPanel(new BorderLayout());
        panel.add(spaceSummaryLabel, BorderLayout.NORTH);
        panel.add(tabs, BorderLayout.CENTER);
        return panel;
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
                    nonHydrogenAtomCache.clear();
                    reactionStatsCache.clear();
                    spaceStatsCache.clear();
                    reactionHistogramModel.setRows(List.of());
                    randomExamplesTableModel.setRows(List.of());
                    reactionSummaryLabel.setText("Select a reaction node to inspect statistics.");
                    spaceHistogramModel.setRows(List.of());
                    largestReactionsTableModel.setRows(List.of());
                    spaceSummaryLabel.setText("Select the space node to inspect global statistics.");
                    rightTopCards.show(rightTopPanel, CARD_SYNTHONS);
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
        SynthonTreeBuilder.NodeData selectedNodeData = nodeData(node);
        List<SynthonTableModel.Row> rows;
        if (selectedNodeData != null && selectedNodeData.type() == SynthonTreeBuilder.NodeType.ROOT) {
            rows = new ArrayList<>();
        } else {
            rows = collectRowsForNode(node);
        }
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
        updateRightTopCard(selectedNodeData);
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

    private void updateRightTopCard(SynthonTreeBuilder.NodeData nodeData) {
        if (nodeData == null) {
            rightTopCards.show(rightTopPanel, CARD_SYNTHONS);
            return;
        }
        if (nodeData.type() == SynthonTreeBuilder.NodeType.ROOT) {
            updateSpaceCard();
            return;
        }
        if (nodeData.type() != SynthonTreeBuilder.NodeType.REACTION) {
            rightTopCards.show(rightTopPanel, CARD_SYNTHONS);
            return;
        }
        rightTopCards.show(rightTopPanel, CARD_REACTION);
        String reactionId = nodeData.reactionId();
        String cacheKey = reactionStatsCacheKey(reactionId);
        ReactionStats cached = reactionStatsCache.get(cacheKey);
        if (cached != null) {
            applyReactionStats(cached);
            return;
        }

        reactionSummaryLabel.setText(String.format("Computing reaction stats for %s ...", reactionId));
        reactionHistogramModel.setRows(List.of());
        randomExamplesTableModel.setRows(List.of());

        SwingWorker<ReactionStats, Void> worker = new SwingWorker<>() {
            @Override
            protected ReactionStats doInBackground() {
                return computeReactionStats(reactionId);
            }

            @Override
            protected void done() {
                try {
                    ReactionStats stats = get();
                    reactionStatsCache.put(cacheKey, stats);

                    DefaultMutableTreeNode selected = selectedNode();
                    SynthonTreeBuilder.NodeData selectedData = nodeData(selected);
                    if (selectedData == null
                            || selectedData.type() != SynthonTreeBuilder.NodeType.REACTION
                            || !Objects.equals(selectedData.reactionId(), reactionId)
                            || !Objects.equals(reactionStatsCacheKey(reactionId), cacheKey)) {
                        return;
                    }
                    applyReactionStats(stats);
                } catch (Exception ex) {
                    showError("Could not compute reaction statistics", ex);
                }
            }
        };
        worker.execute();
    }

    private void updateSpaceCard() {
        rightTopCards.show(rightTopPanel, CARD_SPACE);
        String cacheKey = spaceStatsCacheKey();
        SpaceStats cached = spaceStatsCache.get(cacheKey);
        if (cached != null) {
            applySpaceStats(cached);
            return;
        }

        spaceSummaryLabel.setText("Computing space-level statistics ...");
        spaceHistogramModel.setRows(List.of());
        largestReactionsTableModel.setRows(List.of());

        SwingWorker<SpaceStats, Void> worker = new SwingWorker<>() {
            @Override
            protected SpaceStats doInBackground() {
                return computeSpaceStats();
            }

            @Override
            protected void done() {
                try {
                    SpaceStats stats = get();
                    spaceStatsCache.put(cacheKey, stats);

                    DefaultMutableTreeNode selected = selectedNode();
                    SynthonTreeBuilder.NodeData selectedData = nodeData(selected);
                    if (selectedData == null
                            || selectedData.type() != SynthonTreeBuilder.NodeType.ROOT
                            || !Objects.equals(spaceStatsCacheKey(), cacheKey)) {
                        return;
                    }
                    applySpaceStats(stats);
                } catch (Exception ex) {
                    showError("Could not compute space statistics", ex);
                }
            }
        };
        worker.execute();
    }

    private String reactionStatsCacheKey(String reactionId) {
        return reactionId + "|" + (downsampledSetsCheckBox.isSelected() ? "downsampled" : "full");
    }

    private String spaceStatsCacheKey() {
        return downsampledSetsCheckBox.isSelected() ? "downsampled" : "full";
    }

    private ReactionStats computeReactionStats(String reactionId) {
        RawSynthonSpace.ReactionData reactionData = currentSpace.getReactions().get(reactionId);
        if (reactionData == null) {
            return new ReactionStats(reactionId, 0, 0, 0L, Map.of(), List.of(), "-", 0, 0, 0, 0);
        }

        Map<Integer, List<RawSynthon>> sets = SynthonTreeBuilder.activeSets(reactionData, downsampledSetsCheckBox.isSelected());
        List<Integer> setIndexes = new ArrayList<>(sets.keySet());
        setIndexes.sort(Comparator.naturalOrder());

        int synthonCount = SynthonTreeBuilder.countSynthons(sets);
        long estimatedCompounds = SynthonTreeBuilder.estimateCompoundCount(sets);

        List<Map<Integer, Long>> perSetSizeCounts = new ArrayList<>();
        List<List<RawSynthon>> orderedSets = new ArrayList<>();
        List<String> setSizeSummary = new ArrayList<>();
        for (Integer setIndex : setIndexes) {
            List<RawSynthon> synthons = sets.get(setIndex);
            if (synthons == null || synthons.isEmpty()) {
                continue;
            }
            orderedSets.add(synthons);
            setSizeSummary.add(String.format("S%d=%d", setIndex, synthons.size()));
            Map<Integer, Long> sizeCounts = new HashMap<>();
            for (RawSynthon synthon : synthons) {
                int atoms = estimateSynthonNonHydrogenAtoms(synthon);
                sizeCounts.merge(atoms, 1L, Long::sum);
            }
            perSetSizeCounts.add(sizeCounts);
        }

        Map<Integer, BigInteger> histogram = convolveSizeHistograms(perSetSizeCounts);
        int minAtoms = minHistogramKey(histogram);
        int maxAtoms = maxHistogramKey(histogram);
        int medianAtoms = medianHistogramKey(histogram);
        int modeAtoms = modeHistogramKey(histogram);
        List<RandomExampleRow> randomExamples = sampleRandomExamples(orderedSets, RANDOM_EXAMPLE_COUNT);

        return new ReactionStats(
                reactionId,
                setIndexes.size(),
                synthonCount,
                estimatedCompounds,
                histogram,
                randomExamples,
                String.join(", ", setSizeSummary),
                minAtoms,
                medianAtoms,
                modeAtoms,
                maxAtoms
        );
    }

    private SpaceStats computeSpaceStats() {
        if (currentSpace == null) {
            return new SpaceStats(0, 0, BigInteger.ZERO, 0, 0, 0, 0, List.of(), List.of());
        }

        Map<String, Integer> bucketReactionCounts = initBucketReactionCounts();
        Map<String, BigInteger> bucketCompoundSums = initBucketCompoundSums();
        List<LargestReactionRow> largest = new ArrayList<>();

        int totalReactions = 0;
        int totalSynthons = 0;
        int zeroCompoundReactions = 0;
        int twoSetReactions = 0;
        int threeSetReactions = 0;
        int otherSetReactions = 0;
        BigInteger totalCompounds = BigInteger.ZERO;

        for (Map.Entry<String, RawSynthonSpace.ReactionData> entry : currentSpace.getReactions().entrySet()) {
            String reactionId = entry.getKey();
            RawSynthonSpace.ReactionData reactionData = entry.getValue();
            if (reactionData == null) {
                continue;
            }
            totalReactions++;

            Map<Integer, List<RawSynthon>> sets = SynthonTreeBuilder.activeSets(reactionData, downsampledSetsCheckBox.isSelected());
            int reactionSynthons = SynthonTreeBuilder.countSynthons(sets);
            totalSynthons += reactionSynthons;

            int nonEmptySetCount = countNonEmptySets(sets);
            if (nonEmptySetCount == 2) {
                twoSetReactions++;
            } else if (nonEmptySetCount == 3) {
                threeSetReactions++;
            } else {
                otherSetReactions++;
            }

            long estimatedCompounds = SynthonTreeBuilder.estimateCompoundCount(sets);
            if (estimatedCompounds == 0L) {
                zeroCompoundReactions++;
            }
            BigInteger estimatedBig = BigInteger.valueOf(estimatedCompounds);
            totalCompounds = totalCompounds.add(estimatedBig);

            String bucket = classifySpaceCompoundBucket(estimatedCompounds);
            bucketReactionCounts.merge(bucket, 1, Integer::sum);
            bucketCompoundSums.merge(bucket, estimatedBig, BigInteger::add);

            largest.add(new LargestReactionRow(
                    reactionId,
                    nonEmptySetCount,
                    reactionSynthons,
                    SynthonTreeBuilder.formatMillions(estimatedCompounds),
                    estimatedCompounds
            ));
        }

        largest.sort(
                Comparator.comparingLong(LargestReactionRow::estimatedCompounds).reversed()
                        .thenComparing(LargestReactionRow::reactionId)
        );
        if (largest.size() > SPACE_LARGEST_REACTIONS_LIMIT) {
            largest = new ArrayList<>(largest.subList(0, SPACE_LARGEST_REACTIONS_LIMIT));
        }

        List<SpaceHistogramRow> histogramRows = toSpaceHistogramRows(bucketReactionCounts, bucketCompoundSums, totalReactions);

        return new SpaceStats(
                totalReactions,
                totalSynthons,
                totalCompounds,
                zeroCompoundReactions,
                twoSetReactions,
                threeSetReactions,
                otherSetReactions,
                histogramRows,
                largest
        );
    }

    private void applySpaceStats(SpaceStats stats) {
        String summary = String.format(
                "<html><b>Space Summary</b> | reactions: %d | synthons: %d | enumeratable compounds: %s<br/>set-types: 2-set=%d, 3-set=%d, other=%d | zero-compound reactions: %d</html>",
                stats.totalReactions(),
                stats.totalSynthons(),
                formatBigInteger(stats.totalCompounds()),
                stats.twoSetReactions(),
                stats.threeSetReactions(),
                stats.otherSetReactions(),
                stats.zeroCompoundReactions()
        );
        spaceSummaryLabel.setText(summary);
        spaceHistogramModel.setRows(stats.histogramRows());
        largestReactionsTableModel.setRows(stats.largestReactions());
        statusLabel.setText(String.format(
                "Space stats: reactions=%d, total compounds=%s",
                stats.totalReactions(),
                formatBigInteger(stats.totalCompounds())
        ));
    }

    private Map<String, Integer> initBucketReactionCounts() {
        Map<String, Integer> map = new HashMap<>();
        for (String label : SPACE_BUCKET_LABELS) {
            map.put(label, 0);
        }
        return map;
    }

    private Map<String, BigInteger> initBucketCompoundSums() {
        Map<String, BigInteger> map = new HashMap<>();
        for (String label : SPACE_BUCKET_LABELS) {
            map.put(label, BigInteger.ZERO);
        }
        return map;
    }

    private List<SpaceHistogramRow> toSpaceHistogramRows(Map<String, Integer> bucketReactionCounts,
                                                         Map<String, BigInteger> bucketCompoundSums,
                                                         int totalReactions) {
        List<SpaceHistogramRow> rows = new ArrayList<>(SPACE_BUCKET_LABELS.length);
        for (String label : SPACE_BUCKET_LABELS) {
            int reactions = bucketReactionCounts.getOrDefault(label, 0);
            BigInteger compounds = bucketCompoundSums.getOrDefault(label, BigInteger.ZERO);
            BigDecimal compoundsBillions = new BigDecimal(compounds)
                    .divide(new BigDecimal("1000000000"), 1, RoundingMode.HALF_UP);
            double share = totalReactions == 0 ? 0.0 : (100.0 * reactions / totalReactions);
            rows.add(new SpaceHistogramRow(
                    label,
                    reactions,
                    compoundsBillions.toPlainString(),
                    percentFormat.format(share)
            ));
        }
        return rows;
    }

    private String classifySpaceCompoundBucket(long compounds) {
        if (compounds < 1_000_000L) {
            return "<1M";
        }
        if (compounds < 10_000_000L) {
            return "1M-10M";
        }
        if (compounds < 100_000_000L) {
            return "10M-100M";
        }
        if (compounds < 1_000_000_000L) {
            return "100M-1B";
        }
        if (compounds < 10_000_000_000L) {
            return "1B-10B";
        }
        if (compounds <= 100_000_000_000L) {
            return "10B-100B";
        }
        return ">100B";
    }

    private int countNonEmptySets(Map<Integer, List<RawSynthon>> sets) {
        int count = 0;
        for (List<RawSynthon> synthons : sets.values()) {
            if (synthons != null && !synthons.isEmpty()) {
                count++;
            }
        }
        return count;
    }

    private String formatBigInteger(BigInteger value) {
        if (value == null) {
            return "0";
        }
        String plain = value.toString();
        StringBuilder out = new StringBuilder();
        int length = plain.length();
        for (int i = 0; i < length; i++) {
            if (i > 0 && (length - i) % 3 == 0) {
                out.append(',');
            }
            out.append(plain.charAt(i));
        }
        return out.toString();
    }

    private void applyReactionStats(ReactionStats stats) {
        String summary = String.format(
                "<html><b>%s</b> | %d sets | %d synthons | %sM compounds<br/>Set sizes: %s<br/>Non-H atoms (estimated): min %d, median %d, mode %d, max %d</html>",
                stats.reactionId(),
                stats.setCount(),
                stats.synthonCount(),
                SynthonTreeBuilder.formatMillions(stats.estimatedCompounds()),
                stats.setSizeSummary().isBlank() ? "-" : stats.setSizeSummary(),
                stats.minAtoms(),
                stats.medianAtoms(),
                stats.modeAtoms(),
                stats.maxAtoms()
        );
        reactionSummaryLabel.setText(summary);

        List<ReactionHistogramRow> histogramRows = toHistogramRows(stats.histogram());
        reactionHistogramModel.setRows(histogramRows);
        randomExamplesTableModel.setRows(stats.randomExamples());
        if (!stats.randomExamples().isEmpty()) {
            randomExamplesTable.setRowSelectionInterval(0, 0);
        }
        statusLabel.setText(String.format(
                "Reaction stats: %s (%d histogram bins, %d random examples)",
                stats.reactionId(),
                histogramRows.size(),
                stats.randomExamples().size()
        ));
    }

    private List<ReactionHistogramRow> toHistogramRows(Map<Integer, BigInteger> histogram) {
        if (histogram.isEmpty()) {
            return List.of();
        }
        BigInteger total = BigInteger.ZERO;
        for (BigInteger count : histogram.values()) {
            total = total.add(count);
        }
        List<Integer> sizes = new ArrayList<>(histogram.keySet());
        sizes.sort(Comparator.naturalOrder());

        List<ReactionHistogramRow> rows = new ArrayList<>(sizes.size());
        for (Integer size : sizes) {
            BigInteger count = histogram.getOrDefault(size, BigInteger.ZERO);
            BigDecimal millions = new BigDecimal(count).divide(new BigDecimal("1000000"), 1, RoundingMode.HALF_UP);
            BigDecimal percent = total.signum() == 0
                    ? BigDecimal.ZERO
                    : new BigDecimal(count).multiply(new BigDecimal("100"))
                    .divide(new BigDecimal(total), 4, RoundingMode.HALF_UP);
            rows.add(new ReactionHistogramRow(
                    size,
                    count.toString(),
                    millions.toPlainString(),
                    percentFormat.format(percent.doubleValue())
            ));
        }
        return rows;
    }

    private Map<Integer, BigInteger> convolveSizeHistograms(List<Map<Integer, Long>> perSetSizeCounts) {
        if (perSetSizeCounts.isEmpty()) {
            return Map.of();
        }
        Map<Integer, BigInteger> distribution = new HashMap<>();
        distribution.put(0, BigInteger.ONE);

        for (Map<Integer, Long> setHistogram : perSetSizeCounts) {
            if (setHistogram.isEmpty()) {
                return Map.of();
            }
            Map<Integer, BigInteger> next = new HashMap<>();
            for (Map.Entry<Integer, BigInteger> left : distribution.entrySet()) {
                for (Map.Entry<Integer, Long> right : setHistogram.entrySet()) {
                    int size = left.getKey() + right.getKey();
                    BigInteger add = left.getValue().multiply(BigInteger.valueOf(right.getValue()));
                    next.merge(size, add, BigInteger::add);
                }
            }
            distribution = next;
        }
        return distribution;
    }

    private int minHistogramKey(Map<Integer, BigInteger> histogram) {
        if (histogram.isEmpty()) {
            return 0;
        }
        return Collections.min(histogram.keySet());
    }

    private int maxHistogramKey(Map<Integer, BigInteger> histogram) {
        if (histogram.isEmpty()) {
            return 0;
        }
        return Collections.max(histogram.keySet());
    }

    private int medianHistogramKey(Map<Integer, BigInteger> histogram) {
        if (histogram.isEmpty()) {
            return 0;
        }
        List<Integer> keys = new ArrayList<>(histogram.keySet());
        keys.sort(Comparator.naturalOrder());
        BigInteger total = BigInteger.ZERO;
        for (Integer key : keys) {
            total = total.add(histogram.getOrDefault(key, BigInteger.ZERO));
        }
        BigInteger target = total.add(BigInteger.ONE).divide(BigInteger.valueOf(2));
        BigInteger cumulative = BigInteger.ZERO;
        for (Integer key : keys) {
            cumulative = cumulative.add(histogram.getOrDefault(key, BigInteger.ZERO));
            if (cumulative.compareTo(target) >= 0) {
                return key;
            }
        }
        return keys.get(keys.size() - 1);
    }

    private int modeHistogramKey(Map<Integer, BigInteger> histogram) {
        if (histogram.isEmpty()) {
            return 0;
        }
        int bestKey = 0;
        BigInteger bestValue = BigInteger.valueOf(-1L);
        for (Map.Entry<Integer, BigInteger> entry : histogram.entrySet()) {
            if (entry.getValue().compareTo(bestValue) > 0
                    || (entry.getValue().equals(bestValue) && entry.getKey() < bestKey)) {
                bestKey = entry.getKey();
                bestValue = entry.getValue();
            }
        }
        return bestKey;
    }

    private List<RandomExampleRow> sampleRandomExamples(List<List<RawSynthon>> orderedSets, int targetCount) {
        if (orderedSets.isEmpty()) {
            return List.of();
        }
        for (List<RawSynthon> set : orderedSets) {
            if (set == null || set.isEmpty()) {
                return List.of();
            }
        }

        Set<String> seen = new LinkedHashSet<>();
        List<RandomExampleRow> out = new ArrayList<>();
        Random random = new Random();
        int attempts = 0;
        while (out.size() < targetCount && attempts < RANDOM_SAMPLE_MAX_ATTEMPTS) {
            attempts++;
            List<StereoMolecule> parts = new ArrayList<>(orderedSets.size());
            List<String> partIds = new ArrayList<>(orderedSets.size());
            boolean valid = true;

            for (List<RawSynthon> set : orderedSets) {
                RawSynthon chosen = set.get(random.nextInt(set.size()));
                partIds.add(chosen.getFragmentId() == null ? "n/a" : chosen.getFragmentId());

                StereoMolecule part = parseIdcode(chosen.getIdcode());
                if (part.getAtoms() == 0) {
                    valid = false;
                    break;
                }
                part.setFragment(true);
                parts.add(part);
            }
            if (!valid || parts.isEmpty()) {
                continue;
            }

            try {
                StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts);
                assembled.ensureHelperArrays(Molecule.cHelperCIP);
                String assembledIdcode = assembled.getIDCode();
                if (assembledIdcode == null || assembledIdcode.isBlank() || !seen.add(assembledIdcode)) {
                    continue;
                }
                out.add(new RandomExampleRow(
                        assembledIdcode,
                        countNonHydrogenAtoms(assembled),
                        String.join(" + ", partIds)
                ));
            } catch (Exception ignored) {
                // skip failed random assembly and continue sampling
            }
        }
        return out;
    }

    private int estimateSynthonNonHydrogenAtoms(RawSynthon synthon) {
        String idcode = synthon == null ? null : synthon.getIdcode();
        if (idcode == null || idcode.isBlank()) {
            return 0;
        }
        Integer cached = nonHydrogenAtomCache.get(idcode);
        if (cached != null) {
            return cached;
        }
        int computed = countNonHydrogenAtoms(parseIdcode(idcode));
        nonHydrogenAtomCache.put(idcode, computed);
        return computed;
    }

    private int countNonHydrogenAtoms(StereoMolecule molecule) {
        int count = 0;
        for (int atom = 0; atom < molecule.getAtoms(); atom++) {
            int atomicNo = molecule.getAtomicNo(atom);
            if (atomicNo > 1 && atomicNo < 92) {
                count++;
            }
        }
        return count;
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

    private void updatePreviewFromRandomExampleSelection() {
        int row = randomExamplesTable.getSelectedRow();
        if (row < 0) {
            return;
        }
        RandomExampleRow selected = randomExamplesTableModel.getRow(row);
        setPreviewMolecule(parseIdcode(selected.idcode()));
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

    private record ReactionStats(String reactionId,
                                 int setCount,
                                 int synthonCount,
                                 long estimatedCompounds,
                                 Map<Integer, BigInteger> histogram,
                                 List<RandomExampleRow> randomExamples,
                                 String setSizeSummary,
                                 int minAtoms,
                                 int medianAtoms,
                                 int modeAtoms,
                                 int maxAtoms) {}

    private record SpaceStats(int totalReactions,
                              int totalSynthons,
                              BigInteger totalCompounds,
                              int zeroCompoundReactions,
                              int twoSetReactions,
                              int threeSetReactions,
                              int otherSetReactions,
                              List<SpaceHistogramRow> histogramRows,
                              List<LargestReactionRow> largestReactions) {}

    private record ReactionHistogramRow(int nonHydrogenAtoms,
                                        String compounds,
                                        String compoundsMillions,
                                        String sharePercent) {}

    private record SpaceHistogramRow(String bucket,
                                     int reactions,
                                     String compoundsBillions,
                                     String sharePercent) {}

    private record LargestReactionRow(String reactionId,
                                      int setCount,
                                      int synthonCount,
                                      String compoundsMillions,
                                      long estimatedCompounds) {}

    private record RandomExampleRow(String idcode,
                                    int nonHydrogenAtoms,
                                    String fragmentIds) {}

    private final class ReactionHistogramTableModel extends AbstractTableModel {

        private final String[] columns = {
                "Non-H Atoms",
                "Compounds",
                "Compounds [M]",
                "Share [%]"
        };
        private final List<ReactionHistogramRow> rows = new ArrayList<>();

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
            return columnIndex == 0 ? Integer.class : String.class;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            ReactionHistogramRow row = rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.nonHydrogenAtoms();
                case 1 -> row.compounds();
                case 2 -> row.compoundsMillions();
                case 3 -> row.sharePercent();
                default -> "";
            };
        }

        public void setRows(List<ReactionHistogramRow> newRows) {
            rows.clear();
            rows.addAll(newRows);
            fireTableDataChanged();
        }
    }

    private final class RandomExamplesTableModel extends AbstractTableModel {

        private final String[] columns = {
                "Random Compound[idcode]",
                "Non-H Atoms",
                "Selected Fragments"
        };
        private final List<RandomExampleRow> rows = new ArrayList<>();

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
                case 1 -> Integer.class;
                default -> String.class;
            };
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            RandomExampleRow row = rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.idcode();
                case 1 -> row.nonHydrogenAtoms();
                case 2 -> row.fragmentIds();
                default -> "";
            };
        }

        public void setRows(List<RandomExampleRow> newRows) {
            rows.clear();
            rows.addAll(newRows);
            fireTableDataChanged();
        }

        public RandomExampleRow getRow(int rowIndex) {
            return rows.get(rowIndex);
        }
    }

    private final class SpaceHistogramTableModel extends AbstractTableModel {

        private final String[] columns = {
                "Compounds Per Reaction",
                "Reactions",
                "Compounds [B]",
                "Share [%]"
        };
        private final List<SpaceHistogramRow> rows = new ArrayList<>();

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
                case 1 -> Integer.class;
                default -> String.class;
            };
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SpaceHistogramRow row = rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.bucket();
                case 1 -> row.reactions();
                case 2 -> row.compoundsBillions();
                case 3 -> row.sharePercent();
                default -> "";
            };
        }

        public void setRows(List<SpaceHistogramRow> newRows) {
            rows.clear();
            rows.addAll(newRows);
            fireTableDataChanged();
        }
    }

    private final class LargestReactionsTableModel extends AbstractTableModel {

        private final String[] columns = {
                "Reaction",
                "Sets",
                "Synthons",
                "Compounds [M]"
        };
        private final List<LargestReactionRow> rows = new ArrayList<>();

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
                case 1, 2 -> Integer.class;
                default -> String.class;
            };
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            LargestReactionRow row = rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.reactionId();
                case 1 -> row.setCount();
                case 2 -> row.synthonCount();
                case 3 -> row.compoundsMillions();
                default -> "";
            };
        }

        public void setRows(List<LargestReactionRow> newRows) {
            rows.clear();
            rows.addAll(newRows);
            fireTableDataChanged();
        }
    }

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
