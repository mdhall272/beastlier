/*
* File CategorySet.java
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
package beastlier.outbreak;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.util.Log;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class CategorySet extends BEASTObject {

    final public Input<String> categoryNameInput = new Input<>("categoryName", "name of the category, used as " +
            "metadata name for the outbreak.", Input.Validate.REQUIRED);
    final public Input<String> categoriesInput = new Input<>("value", "categories encoded as case=value pairs " +
            "separated by commas", Input.Validate.REQUIRED);
    final public Input<Outbreak> outbreakInput = new Input<>("taxa", "contains list of ClinicalCases to map " +
            "categories to", Input.Validate.REQUIRED);

    /**
     * String values of categories in order of cases in outbreak*
     */
    protected String[] caseValues;
    protected String[] values;

    Map<String, Integer> map;

    public void initAndValidate() {
        map = new HashMap<>();
        List<String> labels = outbreakInput.get().asStringList();
        String[] categories = categoriesInput.get().split(",");
        caseValues = new String[labels.size()];

        for (String category : categories) {
            category = category.replaceAll("\\s+", " ");
            String[] strs = category.split("=");
            if (strs.length != 2) {
                throw new IllegalArgumentException("could not parse category: " + category);
            }
            String caseID = normalize(strs[0]);
            int caseNr = labels.indexOf(caseID);
            if (caseNr < 0) {
                throw new IllegalArgumentException("Case (" + caseID + ") is not a known case. Spelling error " +
                        "perhaps?");
            }
            caseValues[caseNr] = normalize(strs[1]);
            map.put(caseID, caseNr);

        }

        // sanity check: did we cover all taxa?
        for (int i = 0; i < labels.size(); i++) {
            if (caseValues[i] == null) {
                throw new IllegalArgumentException("no category specified for " + labels.get(i));
            }
        }

        for (int i = 0; i < labels.size(); i++) {
            Log.info.println(labels.get(i) + " = " + caseValues[i]);
        }

    }

    /**
     * some getters and setters *
     */
    public String getCategoryName() {
        return categoryNameInput.get();
    }


    public String getCategoryName(String caseID){
        return caseValues[map.get(caseID)];
    }

    /**
     * remove start and end spaces
     */
    String normalize(String str) {
        if (str.charAt(0) == ' ') {
            str = str.substring(1);
        }
        if (str.endsWith(" ")) {
            str = str.substring(0, str.length() - 1);
        }
        return str;
    }


}
