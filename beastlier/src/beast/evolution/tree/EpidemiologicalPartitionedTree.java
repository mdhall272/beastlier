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

import beast.core.Citation;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.RandomPartition;
import beast.util.Randomizer;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.Outbreak;

import java.util.Arrays;

@Citation(value=
        "Hall M, Woolhouse M, Rambaut A, (2015) Epidemic reconstruction in a \n" +
                "  phylogenetics framework: transmission trees as partitions of the node set. \n" +
                "  PLoS Computational Biology 11(2): e1003537", year = 2015,
        firstAuthorSurname = "hall",
        DOI="10.1371/10.1371/journal.pcbi.1004613")


/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class EpidemiologicalPartitionedTree extends PartitionedTree {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The set of clinical cases");
    public Input<RealParameter> qInput = new Input<>("q", "The set of q parameters for each clinical case; " +
            "rules of the third type only", null, Input.Validate.OPTIONAL);

    private double zeroTime;
    private RealParameter q;
    private Outbreak outbreak;

    private Double[] infectionTimes;
    private Double[] storedInfectionTimes;

    public void initAndValidate(){

        outbreak = outbreakInput.get();

        for(ClinicalCase aCase : outbreak.getEverInfectedCases()){
            if(aCase.getEndTime() > zeroTime){
                zeroTime = aCase.getEndTime();
            }
        }


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

        infectionTimes = new Double[elementList.size()];
        storedInfectionTimes = new Double[elementList.size()];
        Arrays.fill(infectionTimes, null);
        Arrays.fill(storedInfectionTimes, null);
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
        return outbreak.getCaseByID(elementList.get(node.getPartitionElementNumber()));
    }

    public double getInfectionTime(ClinicalCase aCase){

        if(aCase.wasEverInfected()) {
            int partitionElementNumber = elementList.indexOf(aCase.getID());

            PartitionedTreeNode earliestNode = getEarliestNodeInPartition(partitionElementNumber);

            if(!earliestNode.isPartitionDirty() && (rules==Rules.SECOND_TYPE || !q.isDirty(partitionElementNumber)) &&
                    !(earliestNode.isRoot() && rootBranchLength.isDirty(0))){
                return infectionTimes[partitionElementNumber];
            } else {

                double result;

                if (rules == Rules.SECOND_TYPE) {
                    if (!earliestNode.isRoot()) {
                        result = heightToTime(earliestNode.getParent().getHeight());
                    } else {
                        result = heightToTime(earliestNode.getHeight() + getRootBranchLength());
                    }

                } else {
                    if (!earliestNode.isRoot()) {
                        result = heightToTime(earliestNode.getHeight()
                                + q.getValue(partitionElementNumber) * earliestNode.getLength());
                    } else {
                        result = heightToTime(earliestNode.getHeight()
                                + q.getValue(partitionElementNumber) * getRootBranchLength());
                    }
                }
                infectionTimes[partitionElementNumber] = result;

                return result;
            }
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    public Node getEarliestNodeInPartition(ClinicalCase aCase){
        int index = elementList.indexOf(aCase.getID());
        return getEarliestNodeInPartition(index);
    }

    public double getQ(int index){
        return q.getValue(index);
    }

    public ClinicalCase getInfector(ClinicalCase aCase){
        if(aCase.wasEverInfected()) {
            String ancestorElement = getAncestorPartitionElement(aCase.getID());
            return ancestorElement==null ? null : outbreak.getCaseByID(ancestorElement);
        } else {
            return null;
        }
    }

    private ClinicalCase getCase(String id){
        return outbreak.getCaseByID(id);
    }

    // This forces cases in a random starting partition to have infection times at some point while the infector case
    // was infected. The problem comes when the infection branch is longer than its sibling branch and the infection
    // point is at a lower height than the sibling node. This doesn't invariably cause a problem but it can. Here
    // we randomly assign cases for which this happens a qi that has a higher height.

    private void initialiseQs(){
        for(ClinicalCase aCase : outbreak.getEverInfectedCases()){
            int elementNo = getElementNo(aCase);

            Node mrca = getElementMRCA(elementNo);
            if(!mrca.isRoot()) {
                double mrcaLength = mrca.getLength();
                Node sibling = sibling(mrca);
                double siblingLength = sibling.getLength();
                double difference = mrcaLength - siblingLength;
                if(difference > 0 && mrcaLength*q.getValue(elementNo) < difference){
                    q.setValue(elementNo, (difference + Randomizer.nextDouble()*(mrcaLength-difference))/mrcaLength);
                }
            }
        }
    }

    /**
     * Generates a new tree in which the colours along the branches are
     * indicated by the traits of single-child nodes.
     *
     * Caveat: assumes more than one node exists on tree (i.e. leaf != root)
     *
     * @return Flattened tree.
     */
    public Tree getFlattenedTree() {

        // Create new tree to modify.  Note that copy() doesn't
        // initialise the node array lists, so initArrays() must
        // be called manually.
        Tree flatTree = copy();
        flatTree.initArrays();

        int nextNodeNr = getNodeCount();
        PartitionedTreeNode colourChangeNode;

        for (Node node : getNodesAsArray()) {

            int nodeNum = node.getNr();

            Node newNode = flatTree.getNode(nodeNum);

            newNode.setMetaData(elementLabel,
                    ((PartitionedTreeNode) node).getPartitionElementNumber());

            newNode.metaDataString = String.format("%s=\"%s\"", elementLabel,
                    getElementString(((PartitionedTreeNode) node).getPartitionElementNumber()));

        }

        PartitionedTreeNode newRoot = null;

        for(int elementNo = 0; elementNo<elementList.size(); elementNo++){

            PartitionedTreeNode firstNode = getEarliestNodeInPartition(elementNo);
            int firstNo = firstNode.getNr();
            PartitionedTreeNode firstCopy = (PartitionedTreeNode) flatTree.getNode(firstNo);

            PartitionedTreeNode lastNode = (PartitionedTreeNode)firstNode.getParent();

            if(lastNode != null){
                int lastNo = lastNode.getNr();


                PartitionedTreeNode lastCopy = (PartitionedTreeNode) flatTree.getNode(lastNo);

                double oldLength = firstNode.getLength();

                colourChangeNode = new PartitionedTreeNode();
                colourChangeNode.setNr(nextNodeNr);
                colourChangeNode.setID(String.valueOf(nextNodeNr));
                nextNodeNr++;

                firstCopy.setParent(colourChangeNode);
                colourChangeNode.addChild(firstCopy);

                colourChangeNode.setHeight(firstNode.getHeight() + q.getValue(elementNo)*oldLength);
                colourChangeNode.setMetaData(elementLabel, lastNode.getPartitionElementNumber());
                colourChangeNode.setPartitionElementNumber(lastNode.getPartitionElementNumber());

                colourChangeNode.setParent(lastCopy);
                if (lastCopy.getLeft()==firstCopy)
                    lastCopy.setLeft(colourChangeNode);
                else
                    lastCopy.setRight(colourChangeNode);

            } else {
                //dealing with the root

                colourChangeNode = new PartitionedTreeNode();
                colourChangeNode.setNr(nextNodeNr);
                colourChangeNode.setID(String.valueOf(nextNodeNr));
                nextNodeNr++;

                firstCopy.setParent(colourChangeNode);
                colourChangeNode.addChild(firstCopy);

                colourChangeNode.setHeight(firstNode.getHeight() + q.getValue(elementNo)*rootBranchLength.getValue());
                colourChangeNode.setMetaData(elementLabel, -1);
                colourChangeNode.setPartitionElementNumber(-1);

                newRoot = colourChangeNode;
            }

        }
        //and once again this seems the only way to do it.

        flatTree = new Tree(newRoot);
        flatTree.getInternalNodeCount();
        flatTree.getLeafNodeCount();
        flatTree.initAndValidate();

        return flatTree;
    }

    @Override
    public void store(){
        super.store();
        storedInfectionTimes = infectionTimes.clone();
    }

    @Override
    public void restore(){
        super.restore();
        infectionTimes = storedInfectionTimes;
    }

    @Override
    public boolean somethingIsDirty(){
        return super.somethingIsDirty() || (rules==Rules.THIRD_TYPE && q.somethingIsDirty())
                || rootBranchLength.somethingIsDirty();
    }

    public RealParameter getQ(){
        return q;
    }

    @Override
    public void setEverythingDirty(final boolean isDirty) {
        super.setEverythingDirty(isDirty);
        if(isDirty){
            infectionTimes = new Double[elementList.size()];
            storedInfectionTimes = new Double[elementList.size()];
        }
    }


}
