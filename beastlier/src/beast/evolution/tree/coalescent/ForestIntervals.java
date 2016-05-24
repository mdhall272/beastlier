/*
* File ForestIntervals.java
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
package beast.evolution.tree.coalescent;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.*;
import beast.util.HeapSort;

import java.util.AbstractMap;
import java.util.HashMap;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 *
 * A collection of IntervalLists for the connected components of a partitioned tree (Didelot rules or guided by a
 * tree with Cottam rules)
 */

public class ForestIntervals extends CalculationNode {

    public Input<PartitionedTree> treeInput = new Input<>("tree", "The partitioned tree");
    public Input<Boolean> hasDarkAreasInput = new Input<>("darkAreas", "If some areas of the tree are associated with" +
            " unsampled hosts");

    private AbstractMap<Integer,PartitionIntervals> partitionIntervals;
    private boolean hasDarkAreas;

    public void initAndValidate(){

        partitionIntervals = new HashMap<>();
        hasDarkAreas = hasDarkAreasInput.get();

        if(hasDarkAreas){
            partitionIntervals.put(-1, new PartitionIntervals(treeInput.get(), -1));
        }

        for(int i=0; i<treeInput.get().getNElements(); i++){
            partitionIntervals.put(i, new PartitionIntervals(treeInput.get(), i));
        }


    }

    //at present this should only ever be _called_ if it's changed, so can recalculate
    //todo do better than that

    public PartitionIntervals getIntervals(int elementNo){
        PartitionIntervals pi = partitionIntervals.get(elementNo);

        pi.calculateIntervals();

        return pi;
    }




    public class PartitionIntervals extends TreeIntervals{

        private PartitionedTree tree;
        private final int elementNo;
        private int sampleCount;
        private int storedSampleCount;
        private double totalLength;

        private PartitionIntervals(PartitionedTree tree, int elementNo){
            this.elementNo = elementNo;
            this.tree = tree;
            initAndValidate();
        }

        public int getSampleCount() {
            //the number of coalescences, apparently

            if (requiresRecalculation()){
                sampleCount = tree.countNodesInPartition(elementNo, true, true);
            }
            return sampleCount;
        }

        public double getTotalLength(){
            return totalLength;
        }

        public void store(){
            super.store();
            storedSampleCount = sampleCount;
        }

        public void restore(){
            super.restore();
            sampleCount = storedSampleCount;
        }

        /**
         * Recalculates all the intervals for the given beast.tree.
         */
        @SuppressWarnings("unchecked")
        protected void calculateIntervals() {

            int nodeCount = tree.countNodesInPartition(elementNo, false, true);

            //you need one more for the time of infection

            double[] times = new double[nodeCount + 1];
            int[] childCounts = new int[nodeCount + 1];

            ForestIntervals.collectTimes(tree, times, childCounts, elementNo);

            int[] indices = new int[nodeCount + 1];

            HeapSort.sort(times, indices);

            if (intervals == null || intervals.length != nodeCount + 1) {
                intervals = new double[nodeCount + 1];
                lineageCounts = new int[nodeCount + 1];

                storedIntervals = new double[nodeCount + 1];
                storedLineageCounts = new int[nodeCount + 1];

            }

            // start is the time of the first tip
            double start = times[indices[0]];
            int numLines = 0;
            int nodeNo = 0;
            intervalCount = 0;
            while (nodeNo < nodeCount + 1) {

                int lineagesRemoved = 0;
                int lineagesAdded = 0;

                double finish = times[indices[nodeNo]];
                double next;

                do {
                    final int childIndex = indices[nodeNo];
                    final int childCount = childCounts[childIndex];
                    // don't use nodeNo from here on in do loop
                    nodeNo += 1;
                    if (childCount == 0) {
                        lineagesAdded += 1;
                    } else {
                        lineagesRemoved += (childCount - 1);

                        // no mix of removed lineages when 0 th
                        if (multifurcationLimit == 0.0) {
                            break;
                        }
                    }

                    if (nodeNo < nodeCount) {
                        next = times[indices[nodeNo]];
                    } else break;
                } while (Math.abs(next - finish) <= multifurcationLimit);

                if (lineagesAdded > 0) {

                    if (intervalCount > 0 || ((finish - start) > multifurcationLimit)) {
                        intervals[intervalCount] = finish - start;
                        lineageCounts[intervalCount] = numLines;
                        intervalCount += 1;
                    }

                    start = finish;
                }

                // add sample event
                numLines += lineagesAdded;

                if (lineagesRemoved > 0) {

                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;
                    start = finish;
                }
                // coalescent event
                numLines -= lineagesRemoved;

                if(lineagesAdded == 0 && lineagesRemoved == 0){
                    //only happens at the end;
                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;

                }
            }

            totalLength = 0;

            for(int i=0; i<nodeCount+1; i++){
                totalLength += intervals[i];
            }

            intervalsKnown = true;
        }



    }

    protected static void collectTimes(PartitionedTree tree, double[] times, int[] childCounts, int elementNo) {
        if(tree instanceof EpidemiologicalPartitionedTree && tree.rules == PartitionedTree.Rules.DIDELOT){
            int counter = 0;

            for(Node node : tree.getNodesAsArray()){
                PartitionedTreeNode castNode = (PartitionedTreeNode)node;
                if(castNode.getPartitionElementNumber() == elementNo){
                    times[counter] = node.getHeight();
                    childCounts[counter] = node.isLeaf() ? 0 : 2;
                    counter++;
                } else {
                    PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();
                    if(parent != null && parent.getPartitionElementNumber() == elementNo){
                        times[counter] =
                                ((EpidemiologicalPartitionedTree) tree)
                                        .getInfectionHeightByNr(((PartitionedTreeNode) node)
                                                .getPartitionElementNumber());
                        childCounts[counter] = 0;
                        counter++;
                    }
                }
            }

            times[counter] = ((EpidemiologicalPartitionedTree) tree).getInfectionHeightByNr(elementNo);
            childCounts[counter] = 1;

        } else if(tree instanceof GuidedPartitionedTree){
            ((GuidedPartitionedTree) tree).updatePartitions();

            int counter = 0;
            for(Node node : tree.getNodesAsArray()){
                PartitionedTreeNode castNode = (PartitionedTreeNode)node;
                if(castNode.getPartitionElementNumber() == elementNo){
                    times[counter] = node.getHeight();
                    childCounts[counter] = node.isLeaf() ? 0 : 2;
                    counter++;
                }

                int[] elementHistory = (int[])node.getMetaData(GuidedPartitionedTree.elementBranchString);
                double[] heightHistory = (double[])node.getMetaData(GuidedPartitionedTree.heightBranchString);

                if(elementHistory.length>1){
                    if(elementHistory[elementHistory.length-1] == elementNo){
                        times[counter] = heightHistory[elementHistory.length-2];
                        childCounts[counter] = 0;
                        counter++;
                    }
                }

                if(elementHistory.length>2){
                    for(int portion = elementHistory.length-2; portion>=1; portion--){
                        if(elementHistory[portion]==elementNo){
                            times[counter] = heightHistory[portion-1];
                            childCounts[counter] = 0;
                            counter++;
                        }
                    }
                }
            }

            times[counter] =  ((GuidedPartitionedTree) tree)
                    .thisTreeHeight(((GuidedPartitionedTree)tree).getGuideTree().getInfectionHeightByNr(elementNo));
            childCounts[counter] = 1;

        } else {
            throw new RuntimeException("This type of tree has no within-host model.");
        }
    }




}
