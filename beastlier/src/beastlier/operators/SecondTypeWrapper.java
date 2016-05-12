/*
* File SecondTypeWrapper.java
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
package beastlier.operators;

import beast.core.Input;
import beast.evolution.operators.Exchange;
import beast.evolution.operators.SubtreeSlide;
import beast.evolution.operators.TreeOperator;
import beast.evolution.operators.WilsonBalding;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.evolution.tree.Tree;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class SecondTypeWrapper extends TreeOperator {

    public Input<TreeOperator> innerOperatorInput = new Input<>("inner", "The basic phylogeny operator that this " +
            "operator wraps, performing partition changes after the topology proposal");

    private TreeOperator innerOperator;

    Class[] swapOperators = {Exchange.class};
    Class[] sprOperators = {SubtreeSlide.class, WilsonBalding.class};

    public void initAndValidate() {
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules == PartitionedTree.Rules.SECOND_TYPE)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the second type");
        }
        innerOperator = innerOperatorInput.get();
        if(!Arrays.asList(sprOperators).contains(innerOperator.getClass())
                & !Arrays.asList(swapOperators).contains(innerOperator.getClass())){
            throw new RuntimeException("Incompatible inner operator");
        }

    }

    public double proposal() {
        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        double hr = innerOperator.proposal();
        ArrayList<Node> changedNodes = new ArrayList<>();
        for(Node node : tree.getNodesAsArray()){
            if(node.isDirty() == Tree.IS_FILTHY) {
                changedNodes.add(node);
            }
        }
        //Subtree slide in particular may not change the topology, in which case everything is fine
        if(changedNodes.size()!=0){
            if(Arrays.asList(swapOperators).contains(innerOperator.getClass())){
                //this is a node swap operator
                if(changedNodes.size()!=4){
                    throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                            +" behaving in unexpected ways.");
                }

                ArrayList<Node> swappedNodes = new ArrayList<>();
                for(Node node : changedNodes){
                    if(node.getParent().isDirty() == Tree.IS_FILTHY){
                        swappedNodes.add(node);
                    }
                }

                if(swappedNodes.size()!=2){
                    throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                            +" behaving in unexpected ways.");
                }

                int[] nodeElements = new int[2];
                int[] parentElements = new int[2];

                for(int i=0; i<2; i++){
                    PartitionedTreeNode changedNode = ((PartitionedTreeNode)swappedNodes.get(i));
                    nodeElements[i] = changedNode.getPartitionElementNumber();
                    PartitionedTreeNode changedNodeParent = (PartitionedTreeNode)changedNode.getParent();
                    parentElements[i] = changedNodeParent.getPartitionElementNumber();
                }
                if(nodeElements[0]==parentElements[0] || nodeElements[1]==parentElements[1]){
                    //If this is not true there is nothing to do - the result is already a valid partition
                    for(int i=0; i<2; i++){
                        paintUp(nodeElements[1-i], nodeElements[i], (PartitionedTreeNode)swappedNodes.get(i), null);
                    }
                }

            } else if(Arrays.asList(sprOperators).contains(innerOperator.getClass())){
                if(changedNodes.size()!=4 && changedNodes.size()!=5){
                    //five nodes are involved if the root isn't; four are otherwise
                    throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                            +" behaving in unexpected ways.");
                }

                //this is a node transplantation operator.
                //these moves involve:
                //a transplanted node, which now filthy and has a filthy child and a filthy parent if it is not the
                //root
                //the child of the transplanted node that went with it, which is not filthy but has a filthy parent
                //the new parent of the transplanted node, if it exists. If it does it is now filthy and has a filthy
                //child.
                //the new sibling of the transplanted child, which now filthy and has a filthy parent
                //the old parent of the transplanted node, if it existed. if it did it is now filthy and has a filthy
                //child
                //the old sibling of the transplanted child, which now filthy and either has a filthy parent or has
                //become the root

                ArrayList<Node> forIdentification = new ArrayList<>(changedNodes);

                //find the transplantee
                PartitionedTreeNode transplantedNode = null;
                PartitionedTreeNode transplantedChild = null;
                PartitionedTreeNode oldParent = null;
                PartitionedTreeNode oldSibling = null;
                PartitionedTreeNode newParent = null;
                PartitionedTreeNode newSibling = null;

                for(Node node : forIdentification){
                    Node parent = node.getParent();
                    if(parent == null || parent.isDirty() == Tree.IS_FILTHY){
                        for(Node child: node.getChildren()){
                            if(child.isDirty() == Tree.IS_FILTHY){
                                transplantedNode = (PartitionedTreeNode)child;
                                forIdentification.remove(transplantedNode);
                                break;
                            }
                        }
                    }
                }

                if(transplantedNode == null){
                    //I give up
                    throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                            +" behaving in unexpected ways.");
                }
                for(Node child : transplantedNode.getChildren()){
                    if(child.isDirty() == Tree.IS_FILTHY){
                        newSibling = (PartitionedTreeNode)child;
                    } else {
                        transplantedChild = (PartitionedTreeNode)child;
                    }
                }
                for(Node node : forIdentification){
                    if(node == transplantedNode.getParent()){
                        newParent = (PartitionedTreeNode)node;
                        forIdentification.remove(newParent);
                        break;
                    }
                }
                if(forIdentification.size()==1){
                    oldSibling = (PartitionedTreeNode)forIdentification.get(0);
                    forIdentification.remove(oldSibling);
                } else if (forIdentification.size()==2){
                    if(forIdentification.get(1).getParent() == forIdentification.get(0)){
                        oldSibling = (PartitionedTreeNode)forIdentification.get(1);
                        oldParent = (PartitionedTreeNode)forIdentification.get(0);
                        forIdentification.remove(oldParent);
                        forIdentification.remove(oldSibling);
                    } else if(forIdentification.get(0).getParent() == forIdentification.get(1)){
                        oldSibling = (PartitionedTreeNode)forIdentification.get(0);
                        oldParent = (PartitionedTreeNode)forIdentification.get(1);
                        forIdentification.remove(oldParent);
                        forIdentification.remove(oldSibling);
                    } else {
                        //I give up
                        throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                                +" behaving in unexpected ways.");
                    }
                } else {
                    //I give up
                    throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                            +" behaving in unexpected ways.");
                }

                if(forIdentification.size()!=0){
                    //I give up
                    throw new RuntimeException("Operator "+innerOperator.getClass().toString()
                            +" behaving in unexpected ways.");
                }


                //If the child of the moved node is the earliest node in its element:
                if(transplantedChild.getPartitionElementNumber() != transplantedNode.getPartitionElementNumber()){
                    transplantedNode.setPartitionElementNumber(newSibling.getPartitionElementNumber());
                } else {
                    //Change all paintings up the tree from the old child that used to match the moved node to match
                    //the old child
                    paintUp(transplantedNode.getPartitionElementNumber(), oldSibling.getPartitionElementNumber(),
                            oldSibling, transplantedNode);
                    //Change all paintings up the tree from the moved node that used to match the new child to match
                    //the moved node
                    paintUp(newSibling.getPartitionElementNumber(), transplantedNode.getPartitionElementNumber(),
                            transplantedNode, null);
                }
            } else {
                //I don't know what this is
                throw new RuntimeException("Operator class "+innerOperator.getClass()+" not yet " +
                        "supported");
            }
        }
        treeInput.get().setSomethingIsDirty(true);
        return hr;
    }

    private int[] getParentsArray(Tree tree){
        int[] out = new int[tree.getNodeCount()];
        for(int i=0; i<tree.getNodeCount(); i++){
            if(tree.getNode(i)==tree.getRoot()){
                out[i]=-1;
            } else {
                out[i]=tree.getNode(i).getParent().getNr();
            }
        }
        return out;
    }

    private double[] getHeightsArray(Tree tree){
        double[] out = new double[tree.getNodeCount()];
        for(int i=0; i<tree.getNodeCount(); i++){
            out[i]=tree.getNode(i).getHeight();
        }
        return out;
    }

    private void paintUp(int oldElementNo, int newElementNo, PartitionedTreeNode node, PartitionedTreeNode limit){
        if(node.getParent() == null ){
            return;
        }
        PartitionedTreeNode newParent = (PartitionedTreeNode) node.getParent();
        while(newParent!=null && newParent!=limit && newParent.getPartitionElementNumber() == oldElementNo){
            newParent.setPartitionElementNumber(newElementNo);
            newParent = (PartitionedTreeNode) newParent.getParent();
        }
    }

}
