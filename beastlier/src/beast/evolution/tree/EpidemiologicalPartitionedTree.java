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
import beast.util.RandomPartition;
import beast.util.Randomizer;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.Outbreak;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class EpidemiologicalPartitionedTree extends PartitionedTree {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The set of clinical cases");
    public Input<RealParameter> qInput = new Input<>("q", "The set of q parameters for each clinical case; " +
            "rules of the third type only", null, Input.Validate.OPTIONAL);
    public Input<Double> zeroTimeInput = new Input<>("zeroTime", "The time of the last tip on the forwards timescale",
            0.0, Input.Validate.OPTIONAL);

    private double zeroTime;
    private RealParameter q;
    private Outbreak outbreak;

    public void initAndValidate(){

        zeroTime = zeroTimeInput.get();
        outbreak = outbreakInput.get();

        //todo limits for the qs?

        q = qInput.get();

        super.initAndValidate();

        if(q != null && q.getDimension() != elementList.size()){
            throw new IllegalArgumentException("q has the wrong dimension");
        }

        if((rules == Rules.THIRD_TYPE && q==null) || (rules == Rules.SECOND_TYPE && q!=null)){
            throw new IllegalArgumentException("The q parameters should be present if and only if the rules are of" +
                    " the third type");
        }

        // this is a bit of a hack but should make the run start up from a random tree much more often

        if(rules == Rules.THIRD_TYPE && m_initial.get() instanceof RandomPartition){
            initialiseQs();
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
                    PartitionedTreeNode castNode = (PartitionedTreeNode)node;

                    String id = node.getID();
                    if (id != null) {
                        node.setMetaData(traitSet.getTraitName(), (int)traitSet.getValue(id));
                        castNode.setPartitionElementNumber(outbreak.getCases()
                                .indexOf(outbreak.getCase((int)traitSet.getValue(id))));
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


        elementList = new ArrayList<>(outbreak.getEverInfectedCases());

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

    public double getNodeTime(PartitionedTreeNode node){
        return heightToTime(node.getHeight());
    }

    public ClinicalCase getNodeCase(PartitionedTreeNode node){
        return outbreak.getCase(node.getPartitionElementNumber());
    }

    public double getInfectionTime(ClinicalCase aCase){
        if(aCase.wasEverInfected()) {

            int partitionElementNumber = outbreak.getCases().indexOf(aCase);
            int qIndex = outbreak.getEverInfectedCases().indexOf(aCase);

            if (rules == Rules.SECOND_TYPE) {
                return heightToTime(getEarliestNodeInPartition(partitionElementNumber).getHeight());

            } else {
                PartitionedTreeNode mrca = getElementMRCA(partitionElementNumber);
                if (!mrca.isRoot()) {
                    return heightToTime(mrca.getHeight() + q.getValue(qIndex) * mrca.getLength());
                } else {
                    return heightToTime(mrca.getHeight() + q.getValue(qIndex) * rootBranchLength);
                }
            }
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    public Node getEarliestNodeInPartition(ClinicalCase aCase){
        int index = outbreak.getCaseIndex(aCase);
        return getEarliestNodeInPartition(index);
    }

    public double getQ(int index){
        return q.getValue(index);
    }

    public ClinicalCase getInfector(ClinicalCase aCase){
        PartitionedTreeNode elementMRCA = getElementMRCA(outbreak.getCaseIndex(aCase));
        PartitionedTreeNode parent = (PartitionedTreeNode) elementMRCA.getParent();
        if(parent == null){
            return null;
        } else {
            return getNodeCase(parent);
        }
    }

    // This forces cases in a random starting partition to have infection times at some point while the infector case
    // was infected. The problem comes when the infection branch is longer than its sibling branch and the infection
    // point is at a lower height than the sibling node. This doesn't invariably cause a problem but it can. Here
    // we randomly assign cases for which this happens a qi that has a higher height.

    private void initialiseQs(){
        for(ClinicalCase aCase : outbreak.getCases()){
            int caseNo = outbreak.getCaseIndex(aCase);

            Node mrca = getElementMRCA(caseNo);
            if(!mrca.isRoot()) {
                double mrcaLength = mrca.getLength();
                Node sibling = sibling(mrca);
                double siblingLength = sibling.getLength();
                double difference = mrcaLength - siblingLength;
                if(difference > 0 && mrcaLength*q.getValue(caseNo) < difference){
                    q.setValue(caseNo, (difference + Randomizer.nextDouble()*(mrcaLength-difference))/mrcaLength);
                }
            }
        }
    }


}
