/*
* File TransmissionExchangeA.java
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

import beast.core.Description;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;

import java.util.ArrayList;
import java.util.HashSet;


/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */
@Description("Implements branch exchange operations consistent with node partitions, such that the transmission beast.evolution.tree" +
        "is unchanged")
public class TransmissionExchangeA extends TreeOperator {

    @Override
    public void initAndValidate() {
        if(!(treeInput.get(this) instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        double fLogHastingsRatio = 0;

        try {
            fLogHastingsRatio = operate(tree);
        } catch(Exception e){
            return Double.NEGATIVE_INFINITY;
        }

        return fLogHastingsRatio;
    }

    public double operate(final PartitionedTree tree) throws Exception {

        final int nodeCount = tree.getNodeCount();

        PartitionedTreeNode i = (PartitionedTreeNode)tree.getRoot();

        while (i.isRoot()) {
            i = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(nodeCount));
        }

        ArrayList<PartitionedTreeNode> candidates = getPossibleExchanges(tree, i);

        int candidateCount = candidates.size();

        PartitionedTreeNode j = candidates.get(Randomizer.nextInt(candidates.size()));

        int jFirstCandidateCount = getPossibleExchanges(tree, j).size();

        double HRDenom = (1/((double)candidateCount)) + (1/((double)jFirstCandidateCount));

        final Node iP = i.getParent();
        final Node jP = j.getParent();

        exchangeNodes(i, j, iP, jP);

        if( markCladesInput.get() ) {
            Node iup = iP;
            Node jup = jP;
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

        int reverseCandidatesIfirstCount = getPossibleExchanges(tree, i).size();
        int reverseCandidatesJfirstCount = getPossibleExchanges(tree, j).size();

        double HRNum = (1/(double)reverseCandidatesIfirstCount) + (1/(double)reverseCandidatesJfirstCount);

        if(!tree.isValid()) {
            throw new RuntimeException("TEA isn't working properly");
        }

        return Math.log(HRNum/HRDenom);

    }

    public ArrayList<PartitionedTreeNode> getPossibleExchanges(PartitionedTree tree, PartitionedTreeNode node)
            throws Exception{
        ArrayList<PartitionedTreeNode> out = new ArrayList<PartitionedTreeNode>();
        PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();
        if(parent==null){
            throw new RuntimeException("Can't exchange the root node");
        }
        HashSet<PartitionedTreeNode> possibleParentSwaps = tree.getNodesInSameElement(parent);
        for(PartitionedTreeNode newParent : possibleParentSwaps){
            if(!newParent.isLeaf() && newParent!=parent){
                for(int i=0; i<2; i++){
                    PartitionedTreeNode candidate = (PartitionedTreeNode)newParent.getChild(i);
                    if(candidate != parent
                            && node != newParent
                            && candidate.getHeight()  < parent.getHeight()
                            && node.getHeight() < newParent.getHeight() ) {
                        if(out.contains(candidate) || candidate==node){
                            throw new Exception("Adding a candidate that already exists in the list or" +
                                    " the node itself");
                        }
                        out.add(candidate);
                    }
                }
            }
        }
        return out;
    }


    /* exchange sub-trees whose root are i and j */

    protected void exchangeNodes(Node i, Node j,
                                 Node iP, Node jP) {
        // precondition iP -> i & jP -> j
        replace(iP, i, j);
        replace(jP, j, i);
        // postcondition iP -> j & iP -> i
    }
}