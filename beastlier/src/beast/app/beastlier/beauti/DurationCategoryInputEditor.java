/*
* File DurationCategoryInputEditor.java
*
* Copyright (C) 2016 Matthew Hall mdhall@ic.ac.uk
*
* This file is part of BEASTLIER, a BEAST 2 package.
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this program; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package beast.app.beastlier.beauti;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.Parameter;
import beastlier.outbreak.CategorySet;
import beastlier.outbreak.Outbreak;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import java.util.List;


/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class DurationCategoryInputEditor extends InputEditor.Base {

    DurationCategoryTableModel tableModel;
    CategorySet categories;
    Outbreak outbreak;

    public DurationCategoryInputEditor(BeautiDoc doc) {
        super(doc);
    }


    @Override
    public Class<?> type(){
        return CategorySet.class;
    }

    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption,
                     boolean bAddButtons) {

        categories = (CategorySet)input.get();
        categories.initAndValidate();
        outbreak = categories.outbreakInput.get();

        tableModel = new DurationCategoryTableModel(categories);
        JTable table = new JTable(tableModel);


        Box boxVert = Box.createVerticalBox();

        Box boxHoriz = Box.createHorizontalBox();
        boxHoriz.add(Box.createHorizontalGlue());
        boxVert.add(boxHoriz);
        boxVert.add(new JScrollPane(table));

        add(boxVert);


    }

    class DurationCategoryTableModel extends AbstractTableModel {

        CategorySet categories;

        public DurationCategoryTableModel(CategorySet categories) {
            this.categories = categories;
        }

        @Override
        public int getRowCount() {
            return outbreak.getCases().size();
        }

        @Override
        public int getColumnCount() {
            return 2;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            if (columnIndex<0 || columnIndex>=getRowCount())
                return null;

            switch(columnIndex) {
                case 0:
                    // Case name:
                    return outbreak.getCase(rowIndex);
                case 1:
                    // Type:
                    return categories.getDistribution(outbreak.getCase(rowIndex).getID());
                default:
                    return null;
            }
        }

        @Override
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return columnIndex == 1;
        }

        @Override
        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            String caseId = categories.outbreakInput.get().getCase(rowIndex).getID();
            String catString = categories.categoriesInput.get();
            int startIdx = catString.indexOf(caseId + "=");
            int endIdx = catString.indexOf(",", startIdx);

            String newCatString = catString.substring(0, startIdx);
            newCatString += caseId + "=" + aValue;
            if (endIdx>=0)
                newCatString += catString.substring(endIdx);

            categories.categoriesInput.setValue(newCatString, categories);
            try {
                categories.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error setting duration distribution value.");
            }

            fireTableCellUpdated(rowIndex, columnIndex);
        }

        @Override
        public String getColumnName(int column) {
            switch(column) {
                case 0:
                    return "Name";
                case 1:
                    return "Duration Distribution ID";
                default:
                    return null;
            }
        }
    }

    private static String getParameterString(Parameter.Base param) {

        String str = "";
        for (Object value : (List<Object>) param.valuesInput.get()) {
            if (str.length()>0)
                str += " ";
            str += value;
        }

        return str;
    }




}
