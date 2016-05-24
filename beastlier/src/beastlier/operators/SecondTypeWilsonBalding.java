/*
* File SecondTypeWilsonBalding.java
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

import beast.evolution.operators.WilsonBalding;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class SecondTypeWilsonBalding extends WilsonBalding {

    @Override
    public void initAndValidate() {
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules == PartitionedTree.Rules.COTTAM)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the second type");
        }
    }

    @Override
    public double proposal() {
        PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        double oldMinAge, newMinAge, newRange, oldRange, newAge, hastingsRatio;

        // choose a random node avoiding root
        final int nodeCount = tree.getNodeCount();
        PartitionedTreeNode i;
        do {
            i = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot());
        final PartitionedTreeNode iP = (PartitionedTreeNode)i.getParent();

        // choose another random node to insert i above
        PartitionedTreeNode j;
        PartitionedTreeNode jP;

        // make sure that the target branch <k, j> is above the subtree being moved
        do {
            j = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(nodeCount));
            jP = (PartitionedTreeNode)j.getParent();
        } while ((jP != null && jP.getHeight() <= i.getHeight()) || (i.getNr() == j.getNr()));

        // disallow moves that change the root.
        if (j.isRoot() || iP.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        assert jP != null;  // j != root tested above
        final int pnr = iP.getNr();
        final int jPnr = jP.getNr();
        if ( jPnr == pnr || j.getNr() == pnr || jPnr == i.getNr())
            return Double.NEGATIVE_INFINITY;

        final PartitionedTreeNode CiP = (PartitionedTreeNode)getOtherChild(iP, i);

        final PartitionedTreeNode PiP = (PartitionedTreeNode)iP.getParent();

        newMinAge = Math.max(i.getHeight(), j.getHeight());
        newRange = jP.getHeight() - newMinAge;
        newAge = newMinAge + (Randomizer.nextDouble() * newRange);
        oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
        oldRange = PiP.getHeight() - oldMinAge;
        hastingsRatio = newRange / Math.abs(oldRange);

        if (oldRange == 0 || newRange == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
            return Double.NEGATIVE_INFINITY;
        }

        // disconnect p
        final PartitionedTreeNode iG = (PartitionedTreeNode)iP.getParent();
        replace(iG, iP, CiP);
        //iP's subtree is gone, act like it's not there
        if(i.getPartitionElementNumber() == iP.getPartitionElementNumber()) {
            PartitionedTreeNode currentNode = iG;
            while (currentNode != null &&
                    currentNode.getPartitionElementNumber() == iP.getPartitionElementNumber()) {
                currentNode.setPartitionElementNumber(CiP.getPartitionElementNumber());
                currentNode.setPartitionDirty(true);
                currentNode = (PartitionedTreeNode) currentNode.getParent();
            }
        }
        // re-attach, first child node to p
        replace(iP, CiP, j);
        // then parent node of j to p
        replace(jP, j, iP);

        //now it's back
        if(i.getPartitionElementNumber() == iP.getPartitionElementNumber()) {
            PartitionedTreeNode currentNode = jP;
            while (currentNode != null &&
                    currentNode.getPartitionElementNumber() == j.getPartitionElementNumber()) {
                currentNode.setPartitionElementNumber(iP.getPartitionElementNumber());
                currentNode.setPartitionDirty(true);
                currentNode = (PartitionedTreeNode) currentNode.getParent();
            }
        } else {
            iP.setPartitionElementNumber(j.getPartitionElementNumber());
            iP.setPartitionDirty(true);
        }

        // mark paths to common ancestor as changed
        if( markCladesInput.get() ) {
            Node iup = iG;
            Node jup = iP;
            while (iup != jup) {
                if( iup.getHeight() < jup.getHeight() ) {
                    assert !iup.isRoot();
                    iup = iup.getParent();
                    iup.makeDirty(Tree.IS_FILTHY);
                } else {
                    assert !jup.isRoot();
                    jup = jup.getParent();
                    jup.makeDirty(Tree.IS_FILTHY);
                }
            }
        }

//        }

        iP.setHeight(newAge);

        return Math.log(hastingsRatio);
    }

}
