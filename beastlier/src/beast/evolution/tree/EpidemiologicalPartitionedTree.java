/*
* File EpidemiologicalPartitionedTree.java
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
package beast.evolution.tree;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.Outbreak;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class EpidemiologicalPartitionedTree extends PartitionedTree {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The set of clinical cases");
    public Input<List<RealParameter>> qInput = new Input<>("q", "The set of q parameters for each clinical case; " +
            "rules of the third type only", null, Input.Validate.OPTIONAL);
    public Input<Double> zeroTimeInput = new Input<>("zeroTime", "The time of the last tip on the forwards timescale", 0.0,
            Input.Validate.OPTIONAL);

    private double zeroTime;
    private List<RealParameter> q;
    private Outbreak outbreak;

    public void initAndValidate(){
        super.initAndValidate();

        zeroTime = zeroTimeInput.get();
        outbreak = outbreakInput.get();
        q = qInput.get();

        if((rules == Rules.THIRD_TYPE && q==null) || (rules == Rules.SECOND_TYPE && q!=null)){
            throw new IllegalArgumentException("The q parameters should be present if and only if the rules are of" +
                    " the third type");
        }
    }

    @Override
    protected void processTraits(List<TraitSet> traitList) {
        // Record trait set associated with leaf types.
        for (TraitSet traitSet : traitList) {
            if (traitSet.getTraitName().equals(elementLabel)) {
                elementTraitSet = traitSet;
                break;
            }
        }

        if(elementTraitSet == null) {
            throw new RuntimeException("No TraitSet for partition elements specified");
        }


        for (TraitSet traitSet : traitList) {
            if(traitSet != elementTraitSet) {
                for (Node node : getExternalNodes()) {
                    String id = node.getID();
                    if (id != null) {
                        node.setMetaData(traitSet.getTraitName(), traitSet.getValue(id));
                    }
                }
                if (traitSet.isDateTrait())
                    timeTraitSet = traitSet;
            } else {
                for (Node node : getExternalNodes()) {
                    PartitionedTreeNode caseNode = (PartitionedTreeNode)node;

                    String id = node.getID();
                    if (id != null) {
                        node.setMetaData(traitSet.getTraitName(), traitSet.getValue(id));
                        caseNode.setPartitionElementNumber(outbreak.getCases()
                                .indexOf(outbreak.getCaseByID(traitSet.getStringValue(id))));
                    }
                }
            }
        }
        traitsProcessed = true;

        //  todo when implementing BEAUTi
        //
        // Use explicitly-identified type trait set if available.
        // Seems dumb, but needed for BEAUti as ListInputEditors
        // muck things up...
//        if (elementTraitInput.get() != null)
//            elementTraitSet = elementTraitInput.get();


        elementList = new ArrayList<>(outbreak.getCases());

        System.out.println("Partition element trait with the following elements detected:");
        for (int i = 0; i < elementList.size(); i++) {
            System.out.println(elementList.get(i) + " (" + i + ")");
        }

    }


    private double heightToTime(double height){
        return zeroTime - height;
    }

    private double timeToHeight(double time){
        return zeroTime - time;
    }

    private double getInfectionTime(ClinicalCase aCase){
        int partitionElementNumber = outbreak.getCases().indexOf(aCase);

        if(rules==Rules.SECOND_TYPE){
            PartitionedTreeNode tip = getTipsInElement(partitionElementNumber).iterator().next();
            return heightToTime(getEarliestNodeInPartition(partitionElementNumber).getHeight());

        } else {
            PartitionedTreeNode mrca = getElementMRCA(partitionElementNumber);
            if(!mrca.isRoot()) {
                return heightToTime(mrca.getHeight() + q.get(partitionElementNumber).getValue() * mrca.getLength());
            } else {
                return heightToTime(mrca.getHeight() + q.get(partitionElementNumber).getValue() * rootBranchLength);
            }
        }
    }


}
