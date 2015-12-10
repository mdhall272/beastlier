/*
* File Exchange.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
/*
 * ExchangeOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;


/*
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or less tips!
 */

@Description("Implements branch exchange operations consistent with node partitions, also changing the transmission" +
        "tree")
public class TransmissionExchangeB extends TreeOperator {


    @Override
    public void initAndValidate() throws Exception {
        if(!(treeInput.get(this) instanceof PartitionedTree)){
            throw new Exception("This operator is designed for partitioned trees only");
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


    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     * @param tree
     */
    public double operate(final PartitionedTree tree) {

        final int nodeCount = tree.getNodeCount();

        PartitionedTreeNode i = (PartitionedTreeNode)tree.getRoot();
        PartitionedTreeNode iP = null;
        boolean partitionsMatch = true;

        while (i.isRoot() || partitionsMatch) {
            i = (PartitionedTreeNode) tree.getNode(Randomizer.nextInt(nodeCount));
            iP = (PartitionedTreeNode) i.getParent();
            partitionsMatch = i.isRoot() || i.getPartitionElementNumber() == iP.getPartitionElementNumber();
        }

        ArrayList<PartitionedTreeNode> candidates = getPossibleExchanges(tree, i);

        int candidateCount = candidates.size();

        if(candidateCount==0){
            return Double.NEGATIVE_INFINITY;
        }

        PartitionedTreeNode j = candidates.get(Randomizer.nextInt(candidates.size()));

        int jFirstCandidateCount = getPossibleExchanges(tree, j).size();

        double HRDenom = (1/((double)candidateCount)) + (1/((double)jFirstCandidateCount));

        PartitionedTreeNode jP = (PartitionedTreeNode)j.getParent();

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

        return Math.log(HRNum/HRDenom);
    }

    public ArrayList<PartitionedTreeNode> getPossibleExchanges(PartitionedTree tree, PartitionedTreeNode node){
        ArrayList<PartitionedTreeNode> out = new ArrayList<>();
        PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();
        if(parent==null){
            throw new RuntimeException("Can't exchange the root node");
        }
        if(parent.getPartitionElementNumber() == node.getPartitionElementNumber()){
            throw new RuntimeException("This node is not exchangeable by this operator");
        }
        for(int i=0; i<tree.getNodeCount(); i++){

            PartitionedTreeNode candidate = (PartitionedTreeNode)tree.getNode(i);

            PartitionedTreeNode newParent = (PartitionedTreeNode)candidate.getParent();
            if(newParent!=parent && newParent!=null){
                if(candidate != parent
                        && node != newParent
                        && candidate.getHeight() < parent.getHeight()
                        && node.getHeight()  < newParent.getHeight()
                        && newParent.getPartitionElementNumber()!=candidate.getPartitionElementNumber()){
                    if(out.contains(candidate) || candidate==node){
                        throw new RuntimeException("Adding a candidate that already exists in the list or" +
                                " the node itself");
                    }
                    out.add(candidate);
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
