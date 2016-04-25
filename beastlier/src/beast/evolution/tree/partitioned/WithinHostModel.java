/*
* File WithinHostModel.java
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
package beast.evolution.tree.partitioned;

import beast.core.Input;
import beast.evolution.tree.*;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.Outbreak;

import java.util.HashMap;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public abstract class WithinHostModel extends TreeDistribution {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The collection of clinical cases");

    protected EpidemiologicalPartitionedTree tree;
    protected Outbreak outbreak;

    protected HashMap<ClinicalCase, Treelet> elementsAsTrees;
    protected HashMap<ClinicalCase, Treelet> storedElementsAsTrees;

    @Override
    public void initAndValidate(){
        if(!(treeInput.get() instanceof EpidemiologicalPartitionedTree )){
            throw new IllegalArgumentException("Trees given to the within-host model must have node partitions and" +
                    " an outbreak");
        }

        tree = (EpidemiologicalPartitionedTree) treeInput.get();
        outbreak = outbreakInput.get();

        if(tree.getRules() == PartitionedTree.Rules.SECOND_TYPE){
            throw new IllegalArgumentException("Trees must be partitioned by third-type rules to have a within-" +
                    "host model");
        }

        explodeTree();
        storedElementsAsTrees = new HashMap<>(elementsAsTrees);
    }

    protected class Treelet extends Tree {

        private double zeroHeight;

        protected Treelet(Tree tree, double zeroHeight){
            tree.initAndValidate();
            this.zeroHeight = zeroHeight;
        }

        protected double getZeroHeight(){
            return zeroHeight;
        }

        protected void setZeroHeight(double rootBranchLength){
            this.zeroHeight = zeroHeight;
        }
    }

    protected void explodeTree(){

        for(int i=0; i<outbreak.getStateCount(); i++){
            ClinicalCase aCase = outbreak.getCase(i);
            if(aCase.wasEverInfected() && elementsAsTrees.get(aCase)==null){

                Node elementRoot = tree.getEarliestNodeInPartition(aCase);

                double extraHeight;

                if(elementRoot.isRoot()){
                    extraHeight = tree.getRootBranchLength()
                            * tree.getQ(i);
                } else {
                    extraHeight = elementRoot.getLength()
                            * tree.getQ(i);
                }

                Node newRoot = new Node();

                //todo this really might not work

                newRoot.setHeight(0);

                Tree littleTree = new Tree(newRoot);

                if (!elementRoot.isLeaf()) {
                    for (int j = 0; j < elementRoot.getChildCount(); j++) {
                        copyPartitionToTreelet((PartitionedTreeNode)elementRoot.getChild(i), newRoot, aCase);
                    }
                }

                double minHeight = 0;
                for(int nodeNo = 0; nodeNo<littleTree.getNodeCount(); nodeNo++){
                    Node node = littleTree.getNode(nodeNo);
                    if(node.getHeight()<minHeight) {
                        minHeight = node.getHeight();
                    }
                }

                for(int nodeNo = 0; nodeNo<littleTree.getNodeCount(); nodeNo++){
                    Node node = littleTree.getNode(nodeNo);
                    node.setHeight(node.getHeight() + minHeight);
                }


                Treelet treelet = new Treelet(littleTree,
                        littleTree.getRoot().getHeight() + extraHeight);

                elementsAsTrees.put(aCase, treelet);


            }
        }
    }

    private void copyPartitionToTreelet(PartitionedTreeNode oldNode, Node newParent,
                                        ClinicalCase element){

        if (oldNode.getPartitionElementNumber() == outbreak.getCaseIndex(element)) {
            if (oldNode.isLeaf()) {
                Node newTip = new Node(tree.getTaxonId(oldNode));
                newTip.setParent(newParent);
                newTip.setHeight(newParent.getHeight() - oldNode.getLength());
            } else {
                Node newChild = new Node();
                newParent.addChild(newChild);
                newChild.setHeight(newParent.getHeight() - oldNode.getHeight());
                for (int i = 0; i < oldNode.getChildCount(); i++) {
                    PartitionedTreeNode castChild = (PartitionedTreeNode) oldNode.getChild(i);

                    copyPartitionToTreelet(castChild, newChild, element);
                }
            }
        } else {
            // we need a new tip
            Node transmissionTip = new Node("Transmission_" +
                    outbreak.getCase(oldNode.getPartitionElementNumber()).getID());
            double parentTime = tree.getNodeTime((PartitionedTreeNode)oldNode.getParent());
            double childTime = tree.getInfectionTime(tree.getNodeCase(oldNode));
            newParent.addChild(transmissionTip);
            transmissionTip.setHeight(newParent.getHeight() - (childTime - parentTime) );
        }

    }

}
