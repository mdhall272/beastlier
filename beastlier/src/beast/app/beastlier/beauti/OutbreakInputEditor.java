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
import beast.app.draw.ListInputEditor;
import beast.app.util.Utils;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.GeographicallyLocatedClinicalCase;
import beastlier.outbreak.Outbreak;

import javax.swing.*;
import javax.swing.event.CellEditorListener;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.List;

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
    String m_fileName = "";
    String m_caseText = "";
    String caseIDColumnHeader = "Case_ID";
    String endDateColumnHeader = "End_date";
    String latColumnHeader = "Latitude";
    String longColumnHeader = "Longitude";
    JScrollPane scrollPane;
    boolean m_hasGeography = false;
    Object[][] tableData;
    JTable table;

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

        outbreak = (Outbreak)input.get();
        m_hasGeography = (boolean)outbreak.getInputValue("hasGeography");

        JLabel loadFromFileTxt1 = new JLabel("Load epidemiological data from");

        JComboBox fileType = new JComboBox<>(new String[]{".txt", ".csv"});
        fileType.setMaximumSize(new Dimension(50, 24));

        fileType.addActionListener(e -> {
            @SuppressWarnings("unchecked")
            JComboBox<String> combo = (JComboBox<String>) e.getSource();
            m_textFileDelimiter = combo.getSelectedIndex();
        });

        JLabel loadFromFileTxt2 = new JLabel("file:");

        JButton browseButton = new JButton("Browse");

        JButton loadHelp = new JButton("Help");

        Box boxVert = Box.createVerticalBox();

        Box loadRow = Box.createHorizontalBox();

        JTextField txtFile = new JTextField("file");
        txtFile.setMinimumSize(new Dimension(300, 24));
        txtFile.setMaximumSize(new Dimension(300, 24));

        loadRow.add(loadFromFileTxt1);
        loadRow.add(fileType);
        loadRow.add(loadFromFileTxt2);
        loadRow.add(txtFile);
        loadRow.add(browseButton);

        browseButton.addActionListener(e -> {
            File file = Utils.getLoadFile("Load data from file", new File(Beauti.g_sDir), "Select data file", "dat",
                    "txt","csv");
            if (file != null) {
                txtFile.setText(file.getPath());
            }
        });


        loadRow.add(loadHelp);
        loadHelp.addActionListener(e -> {
            JOptionPane.showMessageDialog(boxVert, helpText);
        });

        loadRow.add(Box.createHorizontalGlue());

        JButton loadButton = new JButton("Load");

        loadButton.addActionListener(e -> {
            String textFileDelimiter = ",";

            switch (m_textFileDelimiter){
                case 0:
                    textFileDelimiter = "\t";
                    break;
                case 1:
                    textFileDelimiter = ",";
            }
            try {
                String fileName = txtFile.getText();
                BufferedReader fin = new BufferedReader(new FileReader(txtFile.getText()));
                //the header
                String header = fin.readLine();
                String[] splitHeader = header.split(textFileDelimiter);
                int idCol = -1;
                int endDateCol = -1;
                int latCol = -1;
                int longCol = -1;

                for(int i=0; i<splitHeader.length; i++){
                    if(splitHeader[i].equals(caseIDColumnHeader)){
                        idCol = i;
                    }
                    if(splitHeader[i].equals(endDateColumnHeader)){
                        endDateCol = i;
                    }
                    if(splitHeader[i].equals(latColumnHeader)){
                        latCol = i;
                    }
                    if(splitHeader[i].equals(longColumnHeader)){
                        longCol = i;
                    }
                }


                if (idCol == -1) {
                    JOptionPane.showMessageDialog(this, "Could not find the required header " + caseIDColumnHeader +
                            " in this file");
                }
                if (endDateCol == -1) {
                    JOptionPane.showMessageDialog(this, "Could not find the required header " + endDateColumnHeader +
                            " in this file");
                }
                if ((latCol == -1 & longCol != -1) || (longCol == -1 & latCol != -1)) {
                    JOptionPane.showMessageDialog(this, "Didn't find both latitude and longitude columns");
                }

                List<ClinicalCase> clinicalCaseList = new ArrayList<>();

                while (fin.ready()) {
                    String line = fin.readLine();
                    String[] splitLine = line.split(textFileDelimiter);
                    ClinicalCase aCase;
                    if(latCol!=-1 && longCol!=-1) {
                        aCase = new GeographicallyLocatedClinicalCase();
                        aCase.setID(splitLine[idCol]);
                        if(!(splitLine[endDateCol].equals("NA"))){
                            RealParameter endParam = new RealParameter();
                            endParam.setInputValue("value", Double.parseDouble(splitLine[endDateCol]));
                            aCase.setInputValue("wasEverInfected", true);
                            aCase.setInputValue("endOfInfectiousTime", endParam);
                        } else {
                            aCase.setInputValue("wasEverInfected", false);
                        }
                        RealParameter latParam = new RealParameter();
                        latParam.setInputValue("value", Double.parseDouble(splitLine[latCol]));
                        RealParameter longParam = new RealParameter();
                        longParam.setInputValue("value", Double.parseDouble(splitLine[longCol]));
                        aCase.setInputValue("latitude", latParam);
                        aCase.setInputValue("longitude", longParam);
                    } else {
                        aCase = new ClinicalCase();
                        aCase.setID(splitLine[idCol]);

                        if(!(splitLine[endDateCol].equals("NA"))){
                            RealParameter endParam = new RealParameter();
                            endParam.setInputValue("value", Double.parseDouble(splitLine[endDateCol]));
                            aCase.setInputValue("wasEverInfected", true);
                            aCase.setInputValue("endOfInfectiousTime", endParam);
                        } else {
                            aCase.setInputValue("wasEverInfected", false);
                        }
                    }
                    clinicalCaseList.add(aCase);
                }
                fin.close();
                outbreak.setInputValue("clinicalCase", clinicalCaseList);
                outbreak.setInputValue("hasGeography", latCol!=-1 && longCol!=-1);
                refreshPanel();

            } catch (Exception e2) {
                JOptionPane.showMessageDialog(this, "Loading data from file failed:" + e2.getMessage());
            }
        });

        loadRow.add(loadButton);

        JLabel inputString = new JLabel("Input all case IDs, comma separated:");
        JTextArea inputField = new JTextArea();
        inputField.setMaximumSize(new Dimension(1024, 100));

        JCheckBox geographyCheck =  new JCheckBox("Cases have spatial (lat/long) coordinates");
        geographyCheck.addActionListener(e -> {
            m_hasGeography = geographyCheck.isEnabled();
        });

        JButton goButton = new JButton("Go");

        goButton.addActionListener(e -> {
            String textDelimiter = ",";

            try {
                outbreak = new Outbreak();
                String ids = inputField.getText();
                String[] splitIds = ids.split(textDelimiter);

                ArrayList<ClinicalCase> clinicalCaseList = new ArrayList<ClinicalCase>();

                for(String id : splitIds) {
                    ClinicalCase aCase;
                    if(m_hasGeography) {
                        aCase = new GeographicallyLocatedClinicalCase();
                        aCase.setID(id);
                    } else {
                        aCase = new ClinicalCase();
                        aCase.setID(id);;
                    }
                    aCase.setInputValue("wasEverInfected", true);
                }
                outbreak.setInputValue("clinicalCase", clinicalCaseList);
                outbreak.setInputValue("hasGeography", m_hasGeography);
                refreshPanel();

            } catch (Exception e2) {
                JOptionPane.showMessageDialog(this, "Creation of outbreak element failed: " + e2.getMessage());
            }
        });

        Box stringRow = Box.createHorizontalBox();
        stringRow.add(inputString);
        stringRow.add(inputField);
        stringRow.add(geographyCheck);
        stringRow.add(goButton);

        boxVert.add(loadRow);
        boxVert.add(stringRow);

        if(outbreak.getCases().size()>0){
            boxVert.add(createListBox());
        }

        add(boxVert);
    }

    private Component createListBox() {

        String[] columnData;

        if(m_hasGeography){
            columnData = new String[]{"Name", "Ever infected", "End of infectiousness", "Longitude",
                    "Latitude"};
            tableData = new Object[outbreak.getCases().size()][5];
            convertOutbreakToTableData();
        } else {
            columnData = new String[]{"Name", "Ever infected", "End of infectiousness"};
            tableData = new Object[outbreak.getCases().size()][3];
            convertOutbreakToTableData();

        }

        table = new JTable(tableData, columnData) {
            private static final long serialVersionUID = 1L;

            // method that induces table row shading
            @Override
            public Component prepareRenderer(TableCellRenderer renderer, int Index_row, int Index_col) {
                Component comp = super.prepareRenderer(renderer, Index_row, Index_col);
                //even index, selected or not selected
                if (isCellSelected(Index_row, Index_col)) {
                    comp.setBackground(Color.lightGray);
                } else if (Index_row % 2 == 0 && !isCellSelected(Index_row, Index_col)) {
                    comp.setBackground(new Color(237, 243, 255));
                } else {
                    comp.setBackground(Color.white);
                }
                return comp;
            }
        };

        table.setDefaultEditor(Object.class, new TableCellEditor() {
            JTextField m_textField = new JTextField();
            int m_iRow,
                    m_iCol;

            @Override
            public boolean stopCellEditing() {
                table.removeEditor();
                String text = m_textField.getText();
                tableData[m_iRow][m_iCol] = text;
                convertTableDataToOutbreak();
                convertOutbreakToTableData();
                return true;
            }

            @Override
            public boolean isCellEditable(EventObject anEvent) {
                return table.getSelectedColumn() != 0
                        && (boolean)tableData[table.getSelectedRow()][1];
            }

            @Override
            public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int rowNr,
                                                         int colNr) {
                if (!isSelected) {
                    return null;
                }
                m_iRow = rowNr;
                m_iCol = colNr;
                m_textField.setText((String) value);
                return m_textField;
            }

            @Override
            public boolean shouldSelectCell(EventObject anEvent) {
                return false;
            }

            @Override
            public void removeCellEditorListener(CellEditorListener l) {
            }

            @Override
            public Object getCellEditorValue() {
                return null;
            }

            @Override
            public void cancelCellEditing() {
            }

            @Override
            public void addCellEditorListener(CellEditorListener l) {
            }
        });

        TableColumn col = table.getColumnModel().getColumn(1);
        JCheckBox checkBox = new JCheckBox();
        checkBox.addActionListener(e -> {
            try {
                int row = table.getSelectedRow();
                tableData[row][1] = checkBox.isSelected();
            } catch (Exception ex) {
                ex.printStackTrace();
            }

        });
        col.setCellEditor(new DefaultCellEditor(checkBox));
        col.setCellRenderer(new MyCheckBoxRenderer());
        col.setPreferredWidth(20);
        col.setMaxWidth(20);


        int fontsize = table.getFont().getSize();
        table.setRowHeight(24 * fontsize / 13);
        scrollPane = new JScrollPane(table);

        return scrollPane;
    } // createListBox

    /* synchronise table with data from traitSet BEASTObject */
    private void convertOutbreakToTableData() {
        if(m_hasGeography) {
            for (int i = 0; i < tableData.length; i++) {
                tableData[i][0] = outbreak.getCases().get(i).getID();
                tableData[i][1] = outbreak.getCases().get(i).getInputValue("wasEverInfected");
                if((boolean)outbreak.getCases().get(i).getInputValue("wasEverInfected")) {
                    RealParameter eoitParameter = (RealParameter)outbreak.getCases().get(i)
                            .getInputValue("endOfInfectiousTime");
                    tableData[i][2] = ((List)eoitParameter.getInputValue("value")).get(0);
                } else {
                    tableData[i][2] = "NA";
                }
                RealParameter latParameter = (RealParameter)outbreak.getCases().get(i)
                        .getInputValue("latitude");
                tableData[i][3] = ((List)latParameter.getInputValue("value")).get(0);

                RealParameter longParameter = (RealParameter)outbreak.getCases().get(i)
                        .getInputValue("longitude");
                tableData[i][4] = ((List)longParameter.getInputValue("value")).get(0);

            }
        } else {
            for (int i = 0; i < tableData.length; i++) {
                tableData[i][0] = outbreak.getCases().get(i).getID();
                tableData[i][1] = outbreak.getCases().get(i).getInputValue("wasEverInfected");
                if((boolean)outbreak.getCases().get(i).getInputValue("wasEverInfected")) {
                    RealParameter eoitParameter = (RealParameter)outbreak.getCases().get(i)
                            .getInputValue("endOfInfectiousTime");
                    tableData[i][2] = ((List)eoitParameter.getInputValue("value")).get(0);
                } else {
                    tableData[i][2] = "NA";
                }
            }
        }

        if (table != null) {
            for (int i = 0; i < tableData.length; i++) {
                table.setValueAt(tableData[i][2], i, 2);
                table.setValueAt(tableData[i][3], i, 3);
                table.setValueAt(tableData[i][4], i, 4);
            }
        }
    } // convertOutbreakToTableData

    /**
     * synchronise outbreak BEAST object with table data
     */
    private void convertTableDataToOutbreak() {
        try {
            for (int i = 0; i < tableData.length; i++) {
                outbreak.getCase(i).setInputValue("wasEverInfected", tableData[i][1]);
                if((boolean)tableData[i][1]) {
                    outbreak.getCase(i).setInputValue("endOfInfectiousDate", tableData[i][2]);
                }
                outbreak.getCase(i).setInputValue("longitude", tableData[i][3]);
                outbreak.getCase(i).setInputValue("latitude", tableData[i][4]);
            }
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    public class MyCheckBoxRenderer extends JCheckBox implements TableCellRenderer {
        private static final long serialVersionUID = 1L;

        public MyCheckBoxRenderer() {
            super();
            setOpaque(true);
        }

        @Override
        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column) {
            setSelected((Boolean) value);
            return this;
        }
    }

}
