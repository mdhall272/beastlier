/*
* File TransmissionWilsonBaldingB.java
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
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.util.Randomizer;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class ThirdTypeWilsonBaldingB extends TreeOperator {

    @Override
    public void initAndValidate() {
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules== PartitionedTree.Rules.THIRD_TYPE)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the third type");
        }
    }

    public double proposal() {

        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        PartitionedTreeNode i;

        double oldMinAge, newMinAge, newRange, oldRange, newAge, q;
        // choose a random eligible node
        final int nodeCount = tree.getNodeCount();
        do {
            i = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(nodeCount));
        } while (!eligibleForMove(i, tree));
        final PartitionedTreeNode iP = (PartitionedTreeNode)i.getParent();

        //this one can go anywhere

        PartitionedTreeNode j = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(tree.getNodeCount()));
        PartitionedTreeNode jP = (PartitionedTreeNode)j.getParent();

        while ((jP != null &&  jP.getHeight() <= i.getHeight()) || (i == j)) {
            j = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(tree.getNodeCount()));
            jP = (PartitionedTreeNode)j.getParent();
        }

        if (iP == tree.getRoot() || j == tree.getRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        if (jP == iP || j == iP || jP == i) return Double.NEGATIVE_INFINITY;

        final PartitionedTreeNode CiP = (PartitionedTreeNode)getOtherChild(iP, i);
        PartitionedTreeNode PiP = (PartitionedTreeNode)iP.getParent();

        newMinAge = Math.max(i.getHeight(), j.getHeight());
        newRange = jP.getHeight() - newMinAge;
        newAge = newMinAge + (Randomizer.nextDouble() * newRange);
        oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
        oldRange = PiP.getHeight() - oldMinAge;
        q = newRange / Math.abs(oldRange);

        // need to account for the random repainting of iP

        if(PiP.getPartitionElementNumber() != CiP.getPartitionElementNumber()){
            q *= 0.5;
        }

        if(jP.getPartitionElementNumber() != j.getPartitionElementNumber()){
            q *= 2;
        }
        replace(PiP, iP, CiP);
        // re-attach, first child node to p
        replace(iP, CiP, j);
        // then parent node of j to p
        replace(jP, j, iP);

        iP.setHeight(newAge);

        // repaint the parent to match either its new parent or its new child (50% chance of each).

        if(Randomizer.nextInt(2)==0){
            iP.setPartitionElementNumber(jP.getPartitionElementNumber());
            iP.setMetaData(tree.getElementLabel(), tree.getElementString(jP.getPartitionElementNumber()));
        } else {
            iP.setPartitionElementNumber(j.getPartitionElementNumber());
            iP.setMetaData(tree.getElementLabel(), tree.getElementString(j.getPartitionElementNumber()));
        }

        return Math.log(q);
    }

    private boolean eligibleForMove(PartitionedTreeNode node, PartitionedTree tree){
        // to be eligible for this move, the node's parent must exist and be in a different partition to itself. This
        // forces the parent to be in the same partition as either its grandchild or its other child.

        return (!node.isRoot() && ((PartitionedTreeNode)node.getParent()).getPartitionElementNumber()
                != node.getPartitionElementNumber());
    }

}
