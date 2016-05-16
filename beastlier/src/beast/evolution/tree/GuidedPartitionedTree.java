/*
* File GuidedPartitionedTree.java
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 *
 * A partitioned tree class intended for individual locus phylogenies in a multilocus model. The phylogeny conforms
 * to a transmission structure given by a separate tree with second-type rules, which is the transmission tree.
 */

public class GuidedPartitionedTree extends PartitionedTree {

    public Input<EpidemiologicalPartitionedTree> ttInput = new Input<>("tt", "The transmission tree");

    private EpidemiologicalPartitionedTree tt;

    //what must be added to the height in this tree to get the height in the guide tree
    private double heightAdjustment;

    public void initAndValidate(){
        super.initAndValidate();

        tt = ttInput.get();
        if(tt.getRules() != Rules.SECOND_TYPE){
            throw new IllegalArgumentException("The guide tree must have second-type rules");
        }

        double latestGuideDate = Double.NEGATIVE_INFINITY;

        for(Node node : tt.getExternalNodes()){
            if(tt.getDate(node.getHeight())>latestGuideDate){
                latestGuideDate = tt.getDate(node.getHeight());
            }
        }

        double latestDateHere = Double.NEGATIVE_INFINITY;

        for(Node node : getExternalNodes()){
            if(getDate(node.height)>latestDateHere){
                latestDateHere = getDate(node.height);
            }
        }

        //This is a bad thing that shouldn't happen
        if(latestDateHere > latestGuideDate){
            throw new RuntimeException("Tips in the guide tree cannot be dated earlier than the tips from this tree");
        }

        heightAdjustment = latestGuideDate - latestDateHere;
    }

    private double guideTreeHeight(double thisTreeHeight){
        return thisTreeHeight + heightAdjustment;
    }

    private double thisTreeHeight(double guideTreeHeight){
        return guideTreeHeight - heightAdjustment;
    }

    //Returns the element that the ancestor of this node was present in at the given height (in either this tree or the
    //guide tree.

    private int elementAtHeight(PartitionedTreeNode node, double height, boolean inGuideTree){

        double gtHeight = inGuideTree ? height : guideTreeHeight(height);
        double gtNodeHeight = guideTreeHeight(node.getHeight());

        if(gtHeight < gtNodeHeight){
            throw new IllegalArgumentException("This function is ambiguous if the given height is lower than the " +
                    "node height");
        }

        int nodeElement = node.getPartitionElementNumber();
        PartitionedTreeNode tipInGuideTree = tt.getTip(nodeElement);

        //in the guide tree, internal nodes are transmissions

        PartitionedTreeNode currentNode = tipInGuideTree;
        while(currentNode != null & currentNode.getHeight() < gtHeight){
            currentNode = (PartitionedTreeNode) currentNode.getParent();
        }

        //If the node immediately after this time point is A infecting B, we want to return A
        if(currentNode!=null){
            return currentNode.getPartitionElementNumber();
        }

        //if this height is earlier than the root of the guide tree, need to check the length of its root branch
        //todo At some point this should be a bone fide node, but not today

        if(gtHeight > tt.getRoot().getHeight()+tt.getRootBranchLength() ){
            return -1;
        } else {
            return ((PartitionedTreeNode)tt.getRoot()).getPartitionElementNumber();
        }

    }

    //Partition the internal nodes according to the guide; return false if you can't do it

    private boolean updatePartitions(){

        //this is going to be slow if done like this, but it will catch problems.

        for(Node node : getInternalNodes()){
            PartitionedTreeNode castNode = (PartitionedTreeNode)node;
            List<Integer> requiredElements = new ArrayList<>();

            for(Node tip : node.getAllLeafNodes()){
                PartitionedTreeNode castTip = (PartitionedTreeNode)tip;
                int elementAtNodeHeight = elementAtHeight(castTip, node.getHeight(), false);
                requiredElements.add(elementAtNodeHeight);
            }
            if(requiredElements.size()!=1){
                return false;
            }
            castNode.setPartitionElementNumber(requiredElements.get(0));

        }

        return false;

    }

    public boolean somethingIsDirty(){
        return super.somethingIsDirty() || tt.somethingIsDirty();
    }

    //extracts the subtree taking place in and element and all its ancestors in the transmission tree. Probably
    //extremely slow.

    public PartitionedTree extractAncestralSubtree(int elementNo, HashMap<Node, Node> references){

        List<Integer> ancestralChain = getAncestralChain(elementNo);

        PartitionedTreeNode root = (PartitionedTreeNode)getRoot();

        return new PartitionedTree(copyDown(root, references, ancestralChain));
    }

    private PartitionedTreeNode copyDown(PartitionedTreeNode oldNode, HashMap<Node, Node> references,
                                         List<Integer> chain){
        if(chain.contains(oldNode.getPartitionElementNumber())){
            PartitionedTreeNode copyNode = new PartitionedTreeNode();
            copyNode.height = oldNode.getHeight();
            copyNode.labelNr = oldNode.getNr();
            copyNode.metaDataString = oldNode.metaDataString;
            copyNode.setParent(null);
            copyNode.setID(oldNode.getID());
            copyNode.setPartitionElementNumber(oldNode.getPartitionElementNumber());
            references.put(copyNode, oldNode);

            if(!oldNode.isLeaf()){
                for(int i=0; i<oldNode.getChildCount(); i++){
                    copyNode.addChild(copyDown((PartitionedTreeNode)oldNode.getChild(i), references, chain));
                }
            }

            return copyNode;
        } else {
            //put in a new tip to represent the point at which an infection went outside the chain
            int infectedElement = oldNode.getPartitionElementNumber();
            double newHeight = thisTreeHeight(tt.getInfectionHeightByNr(infectedElement));

            PartitionedTreeNode transNode = new PartitionedTreeNode();
            transNode.height = newHeight;
            transNode.setParent(null);
            transNode.setPartitionElementNumber(infectedElement);
            return transNode;
        }
    }

}
