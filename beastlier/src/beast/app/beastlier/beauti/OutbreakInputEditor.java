/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package beast.app.beastlier.beauti;

import beast.app.beauti.Beauti;
import beast.app.beauti.BeautiDoc;
import beast.app.beauti.GuessPatternDialog;
import beast.app.draw.InputEditor;
import beast.app.draw.ListInputEditor;
import beast.app.util.Utils;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.Outbreak;
import sun.util.calendar.BaseCalendar;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;

/**
 * BEAUti input editor for BEASTLIER outbreaks.
 */
public class OutbreakInputEditor extends ListInputEditor {

//    OutbreakTableModel tableModel;
    Outbreak outbreak;
    String helpText = "<html>The loaded file can be in tab-separated or comma-separated format. The first row " +
            "must <br>contain the column headers, which must include <font color=red>Host_ID</font>. Optional columns "
            + " <font color=red>End_date</font>, <br><font color=red>Latitude</font>, and <font color=red>Longitude</font>" +
            " will also be read. End dates and coordinates must be decimal <br>numbers. Never-infected susceptibles " +
            "should have NA as end date if that column is <br>present.";
    int m_textFileDelimiter = 0;


    public OutbreakInputEditor(BeautiDoc doc) {
        super(doc);

    }

    @Override
    public Class<?> type() {
        return Outbreak.class;
    }

    @Override
    public Class<?> baseType(){
        return ClinicalCase.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption,
                     boolean bAddButtons) {



        JLabel loadFromFileTxt1 = new JLabel("Load epidemiological data from");

        JComboBox fileType = new JComboBox<>(new String[]{".txt", ".csv"});
        fileType.setMaximumSize(new Dimension(50, 24));

        fileType.addActionListener(e -> {
            @SuppressWarnings("unchecked")
            JComboBox<String> combo = (JComboBox<String>) e.getSource();
            m_textFileDelimiter = combo.getSelectedIndex();
        });

        JLabel loadFromFileTxt2 = new JLabel("file:");

        JButton loadButton = new JButton("Browse");

        JButton loadHelp = new JButton("Help");

        Box boxVert = Box.createVerticalBox();

        Box loadRow = Box.createHorizontalBox();


        loadRow.add(loadFromFileTxt1);
        loadRow.add(fileType);
        loadRow.add(loadFromFileTxt2);
        loadRow.add(loadButton);

        loadButton.addActionListener(e -> {
            File file = Utils.getLoadFile("Load data from file", new File(Beauti.g_sDir), "Select data file", "dat",
                    "txt","csv");
            if (file != null) {
                String delimiter;
                switch(m_textFileDelimiter){
                    case 0:
                        delimiter = "\t";
                        break;
                    case 1:
                        delimiter = ",";
                        break;

                }

//                txtFile.setText(file.getPath());
//                readFromFile.setSelected(true);
//                updateFields();
            }
        });


        loadRow.add(loadHelp);
        loadHelp.addActionListener(e -> {
            JOptionPane.showMessageDialog(boxVert, helpText);
        });

        loadRow.add(Box.createHorizontalGlue());

        JLabel inputString = new JLabel("Or input all case IDs, comma separated:");
        JTextArea inputField = new JTextArea();
        inputField.setMaximumSize(new Dimension(1024, 100));
        JButton goButton = new JButton("Go");

        Box stringRow = Box.createHorizontalBox();
        stringRow.add(inputString);
        stringRow.add(inputField);
        stringRow.add(goButton);

        boxVert.add(loadRow);
        boxVert.add(stringRow);

        add(boxVert);
    }

}
