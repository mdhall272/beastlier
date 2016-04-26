/*
* File TransmissionWilsonBaldingA.java
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
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class TransmissionWilsonBaldingA extends TreeOperator{

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
        // choose a random node avoiding root, and nodes that are ineligible for this move because they have nowhere to
        // go

        ArrayList<PartitionedTreeNode> eligibleNodes = getEligibleNodes(tree);

        i = eligibleNodes.get(Randomizer.nextInt(eligibleNodes.size()));

        int eligibleNodeCount = eligibleNodes.size();

        final PartitionedTreeNode iP = (PartitionedTreeNode)i.getParent();
        HashSet<PartitionedTreeNode> nodesInSameElement = tree.getNodesInSameElement(iP);
        HashSet<PartitionedTreeNode> possibleDestinations = new HashSet<>();
        // we can insert the node above OR BELOW any node in the same partition
        for (PartitionedTreeNode node : nodesInSameElement) {
            possibleDestinations.add(node);
            if (!node.isLeaf()) {
                possibleDestinations.add((PartitionedTreeNode)node.getChild(0));
                possibleDestinations.add((PartitionedTreeNode)node.getChild(1));
            }
        }
        PartitionedTreeNode[] pd = possibleDestinations.toArray(new PartitionedTreeNode[possibleDestinations.size()]);

        PartitionedTreeNode j = pd[Randomizer.nextInt(pd.length)];
        PartitionedTreeNode jP = (PartitionedTreeNode)j.getParent();

        while ((jP != null && (jP.getHeight()  <= i.getHeight())) || (i == j)) {
            j = (pd[Randomizer.nextInt(pd.length)]);
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

        replace(PiP, iP, CiP);
        // re-attach, first child node to p
        replace(iP, CiP, j);
        // then parent node of j to p
        replace(jP, j, iP);

        iP.setHeight(newAge);

        // #todo this may need more than copy-and-paste

        if( markCladesInput.get() ) {
            Node iup = PiP;
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

        double logq = Math.log(q);

        int reverseEligibleNodeCount = getEligibleNodes(tree).size();

        logq += Math.log(eligibleNodeCount/reverseEligibleNodeCount);

        return logq;
    }

    private ArrayList<PartitionedTreeNode> getEligibleNodes(PartitionedTree tree){
        ArrayList<PartitionedTreeNode> out = new ArrayList<>();
        for(int i = 0; i < tree.getNodeCount(); i++){
            PartitionedTreeNode node = (PartitionedTreeNode)tree.getNode(i);

            if(eligibleForMove(node, tree)){
                out.add(node);
            }
        }
        return out;
    }


    private boolean eligibleForMove(PartitionedTreeNode node, PartitionedTree tree){
        // to be eligible for this move, the node's parent and grandparent, or parent and sibling, must be in the
        // same partition element (so removing the parent has no effect on the transmission tree)

        PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();
        PartitionedTreeNode grandparent = parent != null ? (PartitionedTreeNode) parent.getParent() : null;
        PartitionedTreeNode sibling = parent == null ? null : (PartitionedTreeNode) getOtherChild(parent, node);

        return  (!node.isRoot() && ((grandparent != null
                && parent.getPartitionElementNumber() == grandparent.getPartitionElementNumber())
                || parent.getPartitionElementNumber() == sibling.getPartitionElementNumber()));
    }

}
