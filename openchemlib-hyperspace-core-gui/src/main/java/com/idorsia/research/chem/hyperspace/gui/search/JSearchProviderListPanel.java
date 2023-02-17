package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.gui.VerticalFlowLayout;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class JSearchProviderListPanel extends JPanel {


    List<AbstractSearchProvider> searchProviders = new ArrayList<>();


    JScrollPane jsp_main;
    JPanel      jp_main;

    public JSearchProviderListPanel() {
        initGUI();
        reinitList();
    }

    public void addSearchProvider(AbstractSearchProvider sp) {

        sp.addSearchProviderListener(new AbstractSearchProvider.SearchProviderListener() {
            @Override
            public void statusChanged() {
                SwingUtilities.invokeLater(new Runnable() {
                    @Override
                    public void run() {
                        reinitList();
                    }
                });
            }

            @Override
            public void requestViewUpdate() {
                SwingUtilities.invokeLater(new Runnable() {
                    @Override
                    public void run() {
                        reinitList();
                    }
                });
            }
        });

        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                searchProviders.add(sp);
                reinitList();
            }
        });
    }

    private void initGUI() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.jp_main = new JPanel();
        this.jp_main.setLayout(new VerticalFlowLayout(VerticalFlowLayout.LEFT,VerticalFlowLayout.TOP));
        this.jsp_main = new JScrollPane(jp_main);

        this.add(this.jsp_main,BorderLayout.CENTER);

        this.setBorder( BorderFactory.createTitledBorder("Search Providers") );
    }

    /**
     * Updates the views from all search providers.
     */
    public void reinitList() {
        this.jp_main.removeAll();
        for(AbstractSearchProvider asp : this.searchProviders) {
            this.jp_main.add( asp.getSearchConfigurationView() );
        }
        this.revalidate();
    }


}
