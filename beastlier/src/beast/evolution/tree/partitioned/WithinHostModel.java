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

import java.util.Arrays;
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

    protected boolean[] treeletsRequiringExtraction;



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

        treeletsRequiringExtraction = new boolean[tree.getNElements()];
        Arrays.fill(treeletsRequiringExtraction, true);

        elementsAsTrees = new HashMap<>();
        explodeTree();
        storedElementsAsTrees = new HashMap<>(elementsAsTrees);
    }

    protected class Treelet extends Tree {

        private double zeroHeight;

        protected Treelet(Tree tree, double zeroHeight){
            tree.initAndValidate();
            assignFrom(tree);
            this.zeroHeight = zeroHeight;
        }

        protected double getZeroHeight(){
            return zeroHeight;
        }

        protected void setZeroHeight(double rootBranchLength){
            this.zeroHeight = zeroHeight;
        }
    }

    // Turn a partitioned phylogeny into a collection of within-host phylogenies

    protected void explodeTree(){

        for(int i=0; i<tree.getElementList().size(); i++){
            if(treeletsRequiringExtraction[i]){
                String caseName = tree.getElementString(i);

                ClinicalCase aCase = outbreak.getCaseByID(caseName);

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

                newRoot.setHeight(0);

                Tree littleTree = new Tree(newRoot);

                if (!elementRoot.isLeaf()) {
                    for (int j = 0; j < elementRoot.getChildCount(); j++) {
                        copyPartitionToTreelet(littleTree, (PartitionedTreeNode)elementRoot.getChild(j), newRoot, aCase);
                    }
                }

                //only way I've found to get a tree to recognise how many nodes it has!

                littleTree = new Tree(newRoot);
                littleTree.getLeafNodeCount();
                littleTree.getInternalNodeCount();

                double minHeight = 0;
                for(int nodeNo = 0; nodeNo<littleTree.getNodeCount(); nodeNo++){
                    Node node = littleTree.getNode(nodeNo);
                    if(node.getHeight()<minHeight) {
                        minHeight = node.getHeight();
                    }
                }

                for(int nodeNo = 0; nodeNo<littleTree.getNodeCount(); nodeNo++){
                    Node node = littleTree.getNode(nodeNo);
                    node.setHeight(node.getHeight() - minHeight);
                }


                Treelet treelet = new Treelet(littleTree,
                        littleTree.getRoot().getHeight() + extraHeight);

                elementsAsTrees.put(aCase, treelet);
            }
        }
    }

    private void copyPartitionToTreelet(Tree protoTreelet, PartitionedTreeNode oldNode, Node newParent,
                                        ClinicalCase element){

        if (oldNode.getPartitionElementNumber() == tree.getElementNo(element)) {
            if (oldNode.isLeaf()) {
                Node newTip = new Node(tree.getTaxonId(oldNode));
                protoTreelet.addNode(newTip);
                newParent.addChild(newTip);
                newTip.setHeight(newParent.getHeight() - oldNode.getLength());
            } else {
                Node newChild = new Node();
                protoTreelet.addNode(newChild);
                newParent.addChild(newChild);
                newChild.setHeight(newParent.getHeight() - oldNode.getLength());
                for (int i = 0; i < oldNode.getChildCount(); i++) {
                    PartitionedTreeNode castChild = (PartitionedTreeNode) oldNode.getChild(i);
                    copyPartitionToTreelet(protoTreelet, castChild, newChild, element);
                }
            }
        } else {
            // we need a new tip
            Node transmissionTip = new Node("Transmission_" +
                    tree.getElementString(oldNode.getPartitionElementNumber()));
            double parentTime = tree.getNodeTime((PartitionedTreeNode)oldNode.getParent());
            double childTime = tree.getInfectionTime(tree.getNodeCase(oldNode));
            protoTreelet.addNode(transmissionTip);
            newParent.addChild(transmissionTip);
            transmissionTip.setHeight(newParent.getHeight() - (childTime - parentTime) );
        }

    }

}
