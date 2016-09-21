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
import beast.app.draw.ListInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.Parameter;
import beast.evolution.tree.TraitSet;
import beastlier.durations.DurationCategory;

import javax.swing.table.AbstractTableModel;
import java.util.ArrayList;
import java.util.List;


/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class DurationCategoryInputEditor extends ListInputEditor {

    List<DurationCategory> categories;

    public DurationCategoryInputEditor(BeautiDoc doc) {
        super(doc);
    }


    @Override
    public Class<?> baseType(){
        return DurationCategory.class;
    }

    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption,
                     boolean bAddButtons) {

        categories = (ArrayList<DurationCategory>)input.get();


    }

    class DurationCategoryTableModel extends AbstractTableModel {

        List<DurationCategory> categories;

        public DurationCategoryTableModel(List<DurationCategory> categories) {
            this.categories = categories;
        }

        @Override
        public int getRowCount() {
            int count = 0;
            for(DurationCategory category : categories){
                count += category.getCases().size();
            }
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
                    // Taxon name:
                    return typeTraitSet.taxaInput.get().getTaxonId(rowIndex);
                case 1:
                    // Type:
                    return typeTraitSet.getStringValue(rowIndex);
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
            String taxon = taxonSet.getTaxonId(rowIndex);
            String traitString = traitSet.traitsInput.get();
            int startIdx = traitString.indexOf(taxon + "=");
            int endIdx = traitString.indexOf(",", startIdx);

            String newTraitString = traitString.substring(0, startIdx);
            newTraitString += taxon + "=" + (String)aValue;
            if (endIdx>=0)
                newTraitString += traitString.substring(endIdx);

            traitSet.traitsInput.setValue(newTraitString, traitSet);
            try {
                traitSet.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error setting type trait value.");
            }

            fireTableCellUpdated(rowIndex, columnIndex);
        }

        @Override
        public String getColumnName(int column) {
            switch(column) {
                case 0:
                    return "Name";
                case 1:
                    return "CaseID";
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
