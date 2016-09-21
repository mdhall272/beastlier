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
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beastlier.durations.DurationDistribution;
import beastlier.durations.FixedValueDurationDistribution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class CategorySet extends BEASTObject {

    final public Input<String> durationNameInput = new Input<>("durationName", "name of the duration (i.e. latent," +
            "infectious), used as metadata name for the outbreak.", Input.Validate.REQUIRED);
    final public Input<String> categoriesInput = new Input<>("value", "categories encoded as case=value pairs " +
            "separated by commas");
    final public Input<Outbreak> outbreakInput = new Input<>("outbreak", "outbreak of ClinicalCases to map " +
            "categories to", Input.Validate.REQUIRED);
    final public Input<List<DurationDistribution>> distributionsInput =  new Input<>("durationDistribution",
            "List of distributions", new ArrayList<>(), Input.Validate.REQUIRED);

    /**
     * String values of categories in order of cases in outbreak*
     */

    Map<String, DurationDistribution> map;

    public void initAndValidate() {
        List<DurationDistribution> durationDistributions = distributionsInput.get();

        map = new HashMap<>();

        List<ClinicalCase> cases = outbreakInput.get().casesInput.get();

        ArrayList<DurationDistribution> toRemove = new ArrayList<>();

        if(cases.size() > 0) {
            if (categoriesInput.get() == null) {
                // since distributionsInput should be specified, give every case the first distribution

                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < outbreakInput.get().casesInput.get().size(); i++) {
                    if (i > 0) {
                        sb.append(",\n");
                    }
                    sb.append(outbreakInput.get().getCase(i).getID());
                    sb.append("=");
                    sb.append(distributionsInput.get().get(0).getID());
                }

                setInputValue("value", sb.toString());
            }


            String[] categories = categoriesInput.get().split(",");
            for (String category : categories) {
                category = category.replaceAll("\\s+", " ");
                String[] strs = category.split("=");
                if (strs.length != 2) {
                    throw new IllegalArgumentException("could not parse category: " + category);
                }
                String caseID = normalize(strs[0]);
                String distributionID = normalize(strs[1]);
                for (DurationDistribution dd : durationDistributions) {
                    if (distributionID.equals(dd.getID())) {
                        map.put(caseID, dd);
                    }
                }
                if(!map.containsKey(caseID)){
                    FixedValueDurationDistribution newD = new FixedValueDurationDistribution();
                    newD.setID(normalize(distributionID));
                    RealParameter newDRP = new RealParameter();
                    newDRP.setInputValue("value", 1.0);
                    newD.setInputValue("length", newDRP);
                    durationDistributions.add(newD);
                    map.put(caseID, newD);
                }
            }
            // sanity check: did we cover all taxa?


            for(ClinicalCase aCase1 : cases) {
                if (map.get(aCase1.getID()) == null) {
                    throw new IllegalArgumentException("no category specified for " + aCase1.getID());
                }
            }

            for(DurationDistribution dd : durationDistributions){
                if(!map.containsValue(dd)){
                    toRemove.add(dd);
                }
            }


            for (ClinicalCase aCase : cases) {
                Log.info.println(aCase.getID() + " = " + map.get(aCase.getID()));
            }
        }
        // clean up

        if(toRemove.size()>0){
            durationDistributions.removeAll(toRemove);
        }

    }

    /**
     * some getters and setters *
     */
    public String getDistributionName() {
        return durationNameInput.get();
    }

    public DurationDistribution getDistribution(String caseID){
        return map.get(caseID);
    }


    public String getDistributionName(String caseID){
        return getDistribution(caseID).getID();
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
