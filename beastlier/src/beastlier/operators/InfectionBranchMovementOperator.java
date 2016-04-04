/*
* File InfectionBranchMovementOperator.java
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

import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.util.Randomizer;

import java.util.HashMap;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class InfectionBranchMovementOperator extends TreeOperator{

    public void initAndValidate() {
        if(!(treeInput.get(this) instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
    }


    public double proposal() {

        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        // find a case whose infection event we are going to move about
        int caseIndexToAdjust = Randomizer.nextInt(tree.getNElements());
        // if the infection event is the seed of the epidemic, we need to try again
        while(((PartitionedTreeNode)tree.getRoot()).getPartitionElementNumber() == caseIndexToAdjust){
            caseIndexToAdjust = Randomizer.nextInt(tree.getNElements());
        }
        // find the child node of the transmission branch
        PartitionedTreeNode earliestNodeInElement = tree.getEarliestNodeInPartition(caseIndexToAdjust);


        double hr = adjustTree(tree, earliestNodeInElement);

        return hr;
    }

    private double adjustTree(PartitionedTree tree, PartitionedTreeNode node) {
        double out;

        PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();

        //this really shouldn't happen, but just in case

        if(parent == null){
            throw new RuntimeException("Can't adjust an infection branch on the root branch");
        }

        int infectedCase = node.getPartitionElementNumber();
        int infectorCase = parent.getPartitionElementNumber();

        PartitionedTreeNode infectedMRCA = tree.getElementMRCA(infectedCase);

        boolean downIsPossible = node != infectedMRCA;
        boolean upIsPossible = !(tree.isRootBlockedBy(infectedCase, infectorCase)
                && tree.isAncestral(parent));

        if(upIsPossible && downIsPossible){
            out = Randomizer.nextBoolean() ? moveUp(tree, node) : moveDown(tree, node);
        } else if(upIsPossible){
            out = moveUp(tree, node);
        } else if(downIsPossible){
            out = moveDown(tree, node);
        } else {
            return Double.NEGATIVE_INFINITY;
        }

        return out;
    }

    private double moveUp(PartitionedTree tree, PartitionedTreeNode node){

        int infectedCase = node.getPartitionElementNumber();

        PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();

        double hr = 0;

        PartitionedTreeNode sibling = (PartitionedTreeNode)getOtherChild(parent, node);

        int infectorCase = parent.getPartitionElementNumber();

        PartitionedTreeNode infectedMRCA = tree.getElementMRCA(infectedCase);
        PartitionedTreeNode infectorMRCA = tree.getElementMRCA(infectorCase);

        HashMap<PartitionedTreeNode, Integer> changesToMake = new HashMap<>();

        if(tree.isAncestral(parent)){

            PartitionedTreeNode grandparent = (PartitionedTreeNode)parent.getParent();
            if(grandparent!=null && grandparent.getPartitionElementNumber() == infectorCase){
                for(PartitionedTreeNode ancestor : tree.samePartitionElementUpTree(parent)){
                    changesToMake.put(ancestor, infectedCase);
                }
            }

            hr += node == infectedMRCA ? Math.log(0.5) : 0;

        } else {

            if(sibling.getPartitionElementNumber() == infectorCase){
                for(PartitionedTreeNode descendant: tree.samePartitionElementDownTree(sibling)){
                    changesToMake.put(descendant, infectedCase);
                }
                changesToMake.put(sibling, infectedCase);
            }

            hr += node == infectedMRCA ? Math.log(0.5) : 0;
        }
        changesToMake.put(parent, infectedCase);

        for(PartitionedTreeNode changingNode : changesToMake.keySet()){
            changingNode.setPartitionElementNumber(infectedCase);
        }

        //todo remove this once you're happy

        if(!tree.isValid()){
            throw new RuntimeException("Tree does not obey the rules");
        }

        //HR adjustments for reverse moves
        if(tree.isAncestral(parent)){
            hr += sibling == infectorMRCA ? Math.log(2) : 0;
        } else {
            PartitionedTreeNode grandparent = (PartitionedTreeNode)parent.getParent();

            hr += tree.isRootBlockedBy(infectedCase, infectorCase)
                    && tree.isAncestral(grandparent) ? Math.log(2) : 0;

        }

        return hr;
    }

    private double moveDown(PartitionedTree tree, PartitionedTreeNode node){

        int infectedCase = node.getPartitionElementNumber();

        PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();

        int infectorCase = parent.getPartitionElementNumber();

        HashMap<PartitionedTreeNode, Integer> changesToMake = new HashMap<>();

        double out = 0;

        PartitionedTreeNode infectedMRCA = tree.getElementMRCA(infectedCase);

        // check if either child is not ancestral (at most one is not, and if so it must have been in the same
        // partition as both the other child and 'node')
        for(Node child : node.getChildren()){
            PartitionedTreeNode castChild = (PartitionedTreeNode)child;

            if(tree.isAncestral(castChild)){

                for(PartitionedTreeNode descendant: tree.samePartitionElementDownTree(castChild)){
                    changesToMake.put(descendant, infectorCase);
                }
                changesToMake.put(castChild, infectorCase);
            } else if(castChild == infectedMRCA && castChild.getPartitionElementNumber() == infectedCase){
                // we're moving a transmission event as far down as it can go and need to adjust the HR accordingly
                out += Math.log(2);
            }
        }

        //if you couldn't move it any further up
        out += tree.isRootBlockedBy(infectedCase, infectorCase) && tree.isAncestral(parent) ? Math.log(0.5) : 0;

        changesToMake.put(node, infectorCase);

        for(PartitionedTreeNode changingNode : changesToMake.keySet()){
            changingNode.setPartitionElementNumber(infectedCase);
        }

        //todo remove this once you're happy

        if(!tree.isValid()){
            throw new RuntimeException("Tree does not obey the rules");
        }

        return out;
    }


}
