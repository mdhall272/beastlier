/*
* File TransmissionSubtreeSlideB.java
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
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class ThirdTypeSubtreeSlideB extends TreeOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta",
            true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is " +
            "automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by tree-height/log2(n-taxa).", -1.0);
    // shadows size
    double size;
    private double limit;

    @Override
    public void initAndValidate() {
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules== PartitionedTree.Rules.DIDELOT)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the third type");
        }
        size = sizeInput.get();
        limit = limitInput.get();
    }

    public double proposal() {

        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        double logq;

        PartitionedTreeNode i;

        // 1. choose a random eligible node avoiding root
        do {
            i = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(tree.getNodeCount()));
        } while (!eligibleForMove(i, tree));

        final PartitionedTreeNode iP = (PartitionedTreeNode)i.getParent();
        final PartitionedTreeNode CiP = (PartitionedTreeNode)getOtherChild(iP, i);
        final PartitionedTreeNode PiP = (PartitionedTreeNode)iP.getParent();

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = iP.getHeight();
        final double newHeight = oldHeight + delta;

        int iCase = i.getPartitionElementNumber();
        int iPCase = iP.getPartitionElementNumber();
        int CiPCase = CiP.getPartitionElementNumber();
        int PiPCase = -1;
        if(PiP!=null){
            PiPCase = PiP.getPartitionElementNumber();
        }



        // 3. if the move is down
        if (delta > 0) {


            // 3.1 if the topology will change
            if (PiP != null && PiP.getHeight() < newHeight) {
                // find new parent
                PartitionedTreeNode newParent = PiP;
                PartitionedTreeNode newChild = iP;
                while (newParent.getHeight() < newHeight) {
                    newChild = newParent;
                    newParent = (PartitionedTreeNode)newParent.getParent();
                    if (newParent == null) break;
                }

                // 3.1.1 if creating a new root
                if (newChild.isRoot()) {
                    replace(iP, CiP, newChild);
                    replace(PiP, iP, CiP);

                    iP.setParent(null);
                    tree.setRoot(iP);

                }
                // 3.1.2 no new root
                else {
                    replace(iP, CiP, newChild);
                    replace(PiP, iP, CiP);
                    replace(newParent, newChild, iP);
                }

                iP.setHeight(newHeight);


                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChild, oldHeight, null);

                logq = -Math.log(possibleSources);

                // Randomly assign iP the partition of either its parent or the child that is not i, and adjust q
                // appropriately

                if(PiPCase != CiPCase){
                    logq += Math.log(0.5);
                }

                int newiPCase;

                int newChildCase = newChild.getPartitionElementNumber();

                if(newParent != null && newParent.getPartitionElementNumber() != newChildCase){
                    if(Randomizer.nextInt(2)==0){
                        newiPCase = newParent.getPartitionElementNumber();
                    } else {
                        newiPCase = newChildCase;
                    }

                    logq += Math.log(2);
                } else {
                    newiPCase = newChildCase;
                }

                iP.setPartitionElementNumber(newiPCase);
                iP.setMetaData(tree.getElementLabel(), tree.getElementString(newiPCase));

            } else {
                // just change the node height
                // todo you could actually randomise whether the subtree containing iP is changed here

                iP.setHeight(newHeight);
                logq = 0.0;
            }
        }
        // 4 if we are sliding the subtree up.
        else {

            // 4.0 is it a valid move?
            if (i.getHeight() > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology?
            if (CiP.getHeight() > newHeight) {

                List<Node> newChildren = new ArrayList<>();
                final int possibleDestinations = intersectingEdges(CiP, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = Randomizer.nextInt(newChildren.size());
                PartitionedTreeNode newChild = (PartitionedTreeNode)newChildren.get(childIndex);
                PartitionedTreeNode newParent = (PartitionedTreeNode)newChild.getParent();

                // 4.1.1 if iP was root
                if (iP.isRoot()) {
                    // new root is CiP
                    replace(iP, CiP, newChild);
                    replace(newParent, newChild, iP);

                    CiP.setParent(null);
                    tree.setRoot(CiP);

                } else {
                    replace(iP, CiP, newChild);
                    replace(PiP, iP, CiP);
                    replace(newParent, newChild, iP);
                }

                iP.setHeight(newHeight);

                logq = Math.log(possibleDestinations);

                // Randomly assign iP the partition of either its parent or the child that is not i, and adjust q
                // appropriately

                if(PiP!=null && PiPCase != CiPCase){
                    logq += Math.log(0.5);
                }

                int newiPCase;

                int newChildCase = newChild.getPartitionElementNumber();

                if(newParent.getPartitionElementNumber() != newChildCase){
                    if(Randomizer.nextInt(2)==0){
                        newiPCase = newParent.getPartitionElementNumber();
                    } else {
                        newiPCase = newChildCase;
                    }

                    logq += Math.log(2);
                } else {
                    //upward, so don't have to worry about newParent being the root if the topology changed
                    newiPCase = newChildCase;
                }

                iP.setPartitionElementNumber(newiPCase);
                iP.setMetaData(tree.getElementLabel(), tree.getElementString(newiPCase));

            } else {
                iP.setHeight(newHeight);
                logq = 0.0;
            }
        }

        if (logq == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;

        return logq;
    }

    private boolean eligibleForMove(PartitionedTreeNode node, PartitionedTree tree){
        // to be eligible for this move, the node's parent must exist and be in a different partition to itself. This
        // forces the parent to be in the same partition as either its grandchild or its child.

        return (!node.isRoot() && ((PartitionedTreeNode)node.getParent()).getPartitionElementNumber()
                != node.getPartitionElementNumber());
    }

    //intersectingEdges is the same as in normal STS, since there's no additional restriction in this case regarding
    // where nodes can go, and the move does not modify eligibility for itself

    private int intersectingEdges(Node node, double height, List<Node> directChildren) {
        final Node parent = node.getParent();
        if (parent.getHeight() < height) return 0;
        if (node.getHeight() < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }
        if (node.isLeaf()) {
            return 0;
        } else {
            final int count = intersectingEdges(node.getLeft(), height, directChildren) +
                    intersectingEdges(node.getRight(), height, directChildren);
            return count;
        }
    }

    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
//            double f = Math.exp(delta);
            if( limit > 0 ) {
                final Tree tree = treeInput.get();
                final double h = tree.getRoot().getHeight();
                final double k = Math.log(tree.getLeafNodeCount()) / Math.log(2);
                final double lim = (h / k) * limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
                size = f;
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }

}
