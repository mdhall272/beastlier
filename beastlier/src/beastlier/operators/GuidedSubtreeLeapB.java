/*
* File GuidedTransmissionSubtreeLeap.java
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
import beast.evolution.tree.*;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 *
 * This performs the subtree leap move on the transmission tree, and rewires every locus phylogeny accordingly
 */

public class GuidedSubtreeLeapB extends TreeOperator {
    public Input<List<GuidedPartitionedTree>> phylogeniesInput = new Input<>("phylogenies", "All the locus " +
            "phylogenies");
    public Input<Double> sizeInput = new Input<>("size", "");

    private List<GuidedPartitionedTree> phylogenies;
    private double size;

    public void initAndValidate(){
        if(((PartitionedTree)treeInput.get()).getRules() != PartitionedTree.Rules.SECOND_TYPE){
            throw new IllegalArgumentException("The tree must be partitioned according to second-type rules");
        }

        phylogenies = phylogeniesInput.get();
        for(PartitionedTree pTree : phylogenies) {
            if (!(pTree instanceof GuidedPartitionedTree)) {
                throw new IllegalArgumentException("Phylogenies must be guided for this operator");
            }
        }

        size = sizeInput.get();
    }


    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double proposal() {

        EpidemiologicalPartitionedTree tTree = (EpidemiologicalPartitionedTree)treeInput.get();
        double logq;

        final double delta = getDelta();

        final Node root = tTree.getRoot();

        PartitionedTreeNode node;

        do {
            // choose a random node avoiding root
            node = (PartitionedTreeNode)tTree.getNode(Randomizer.nextInt(tTree.getNodeCount()));

        } while (node == root);

        // get its parent - this is the node we will prune/graft
        final PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();

        // get the node's sibling
        final PartitionedTreeNode sibling = (PartitionedTreeNode)getOtherChild(parent, node);

        // and its grand parent
        final PartitionedTreeNode grandParent = (PartitionedTreeNode)parent.getParent();

        final Map<PartitionedTreeNode, Double> destinations = getDestinations(node, parent, sibling, delta);
        final List<PartitionedTreeNode> destinationNodes = new ArrayList<>(destinations.keySet());

        // pick uniformly from this list
        int r = Randomizer.nextInt(destinations.size());

        double forwardProbability = 1.0 / destinations.size();

        final PartitionedTreeNode j = destinationNodes.get(r);
        final double newHeight = destinations.get(j);

        final PartitionedTreeNode jParent = (PartitionedTreeNode)j.getParent();

        if (jParent != null && newHeight > jParent.getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < j.getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        double adjustment = 0;

        if (j == parent || jParent == parent) {
            // the subtree is not actually moving but the height will change

            // todo at present this doesn't affect the locus phylogenies at all and it could cause them to be
            // todo invalid if an infection time is pushed too far down. Can be fixed by moving the hooks with the
            // todo node

            adjustment = 0;
        } else {
            adjustment = unhookRehook(node, j, newHeight);
        }

        parent.setHeight(newHeight);

        if (parent.getParent() != null && newHeight > parent.getParent().getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < node.getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < getOtherChild(parent, node).getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        final Map<PartitionedTreeNode, Double> reverseDestinations = getDestinations(node, parent,
                (PartitionedTreeNode)getOtherChild(parent, node), delta);
        double reverseProbability = 1.0 / reverseDestinations.size();

        // hastings ratio = reverse Prob / forward Prob
        logq = Math.log(reverseProbability) - Math.log(forwardProbability);
        logq += adjustment;
        return logq;
    }

    private Map<PartitionedTreeNode, Double> getDestinations(PartitionedTreeNode node, PartitionedTreeNode parent,
                                                 PartitionedTreeNode sibling, double delta) {

        final Map<PartitionedTreeNode, Double> destinations = new HashMap<PartitionedTreeNode, Double>();

        // get the parent's height
        final double height = parent.getHeight();

        final double heightBelow = height - delta;

        if (heightBelow > node.getHeight()) {
            // the destination height below the parent is compatible with the node
            // see if there are any destinations on the sibling's branch
            final List<PartitionedTreeNode> edges = new ArrayList<PartitionedTreeNode>();

            getIntersectingEdges(sibling, heightBelow, edges);

            // add the intersecting edges and the height
            for (PartitionedTreeNode n : edges) {
                destinations.put(n, heightBelow);
            }
        }

        final double heightAbove = height + delta;

        PartitionedTreeNode node1 = parent;

        // walk up to root
        boolean done = false;
        while (!done) {
            PartitionedTreeNode parent1 = (PartitionedTreeNode)node1.getParent();

            if (parent1 != null) {
                final double height1 = parent1.getHeight();
                if (height1 < heightAbove) {
                    // haven't reached the height above the original height so go down
                    // the sibling subtree
                    PartitionedTreeNode sibling1 = (PartitionedTreeNode)getOtherChild(parent1, node1);

                    double heightBelow1 = height1 - (heightAbove - height1);

                    if (heightBelow1 > node.getHeight()) {

                        final List<PartitionedTreeNode> edges = new ArrayList<PartitionedTreeNode>();

                        getIntersectingEdges(sibling1, heightBelow1, edges);

                        // add the intersecting edges and the height
                        for (PartitionedTreeNode n : edges) {
                            destinations.put(n, heightBelow1);
                        }
                    }
                } else {
                    // add the current node as a destination
                    destinations.put(node1, heightAbove);
                    done = true;
                }

                node1 = parent1;
            } else {
                // node1 is the root - add it as a destination and stop loop
                destinations.put(node1, heightAbove);
                done = true;
            }
        }

        return destinations;
    }

    private double getDelta() {
        return Math.abs(Randomizer.nextGaussian() * size);
    }

    private int getIntersectingEdges(PartitionedTreeNode node, double height,
                                     List<PartitionedTreeNode> edges) {

        final PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();

        if (parent.getHeight() < height) return 0;

        if (node.getHeight() < height) {
            edges.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < node.getChildCount(); i++) {
            count += getIntersectingEdges((PartitionedTreeNode)node.getChild(i), height, edges);
        }
        return count;
    }

    //Do the move and make the respective adjustment to each locus phylogeny. The return value is the adjustment to
    // the HR. This should be called before any attempt to recolour the nodes.

    //For clarity, rather than speed:
    // 1) Remove the moving node in the transmission tree
    // 2) unhook everything
    // 3) calculate reverse probabilities
    // 4) identify destinations and calculate forwards probablities
    // 5) rehook everything
    // 6) reattach the moving node in the transmission tree

    private double unhookRehook(PartitionedTreeNode i, PartitionedTreeNode j, double newTTheight){

        final EpidemiologicalPartitionedTree tt =  (EpidemiologicalPartitionedTree)treeInput.get();

        // get its parent - this is the node we will prune/graft
        final PartitionedTreeNode iP = (PartitionedTreeNode)i.getParent();

        // get the node's sibling
        final PartitionedTreeNode iS = (PartitionedTreeNode)getOtherChild(iP, i);

        // and its grand parent
        final PartitionedTreeNode iG = (PartitionedTreeNode)iP.getParent();

        final PartitionedTreeNode jP = (PartitionedTreeNode)j.getParent();

        double hrAdjust = 0;

        if (iG == null) {
            // if the parent of the original node is the root then the sibling becomes
            // the root.
            iP.removeChild(iS);
            treeInput.get().setRoot(iS);

        } else {
            // remove the parent of node by connecting its sibling to its grandparent.
            iP.removeChild(iS);
            iG.removeChild(iP);
            iG.addChild(iS);
        }

        //the parent's subtree is gone, act like it's not there
        if(i.getPartitionElementNumber() == iP.getPartitionElementNumber()) {
            PartitionedTreeNode currentNode = iG;
            while (currentNode != null &&
                    currentNode.getPartitionElementNumber() == iP.getPartitionElementNumber()) {
                currentNode.setPartitionElementNumber(iS.getPartitionElementNumber());
                currentNode.setPartitionDirty(true);
                currentNode = (PartitionedTreeNode) currentNode.getParent();
            }
        }

        for(GuidedPartitionedTree tree : phylogenies){
            tree.updatePartitions();
        }

        for(GuidedPartitionedTree tree : phylogenies){
            List<Double> heightsToCheck = new ArrayList<>();

            List<PartitionedTreeNode> hooks = tree.getHooks(i.getPartitionElementNumber(),
                    tree.thisTreeHeight(iP.getHeight()));

            for(PartitionedTreeNode hook : hooks){
                PartitionedTreeNode hookParent = (PartitionedTreeNode)hook.getParent();
                if(hookParent != null){
                    heightsToCheck.add(hookParent.getHeight());
                    PartitionedTreeNode hookSibling = (PartitionedTreeNode)getOtherChild(hookParent, hook);
                    PartitionedTreeNode hookGrandparent = (PartitionedTreeNode)hookParent.getParent();
                    if(hookGrandparent != null){
                        hookParent.removeChild(hookSibling);
                        hookGrandparent.removeChild(hookParent);
                        hookGrandparent.addChild(hookSibling);
                    } else {
                        hookParent.removeChild(hookSibling);
                        tree.setRoot(hookSibling);
                    }
                } else {
                    //Actually, I'm reasonably sure this never happens so long as every element is present on every
                    //tree. Assume that for now.
                }
            }

            //this adjusts the HR for the number of places the hook could have been if the move was reversed

            int currentBranchElement = iS.getPartitionElementNumber();
            for(Double height : heightsToCheck){
                while(height > tree.thisTreeHeight(tt.getInfectionHeightByNr(currentBranchElement))){
                    currentBranchElement = tt.getAncestorPartitionElement(currentBranchElement);
                }

                hrAdjust -= Math.log(tree.getBundle((PartitionedTreeNode)tree.getRoot(), height, currentBranchElement,
                        null));
            }

            for(PartitionedTreeNode hook : hooks){
                PartitionedTreeNode hookParent = (PartitionedTreeNode)hook.getParent();
                if(hookParent != null){
                    heightsToCheck.add(hookParent.getHeight());
                    PartitionedTreeNode hookSibling = (PartitionedTreeNode)getOtherChild(hookParent, hook);
                    PartitionedTreeNode hookGrandparent = (PartitionedTreeNode)hookParent.getParent();
                    if(hookGrandparent != null){
                        hookParent.removeChild(hookSibling);
                        hookGrandparent.removeChild(hookParent);
                        hookGrandparent.addChild(hookSibling);
                    } else {
                        hookParent.removeChild(hookSibling);
                        tree.setRoot(hookSibling);
                    }
                } else {
                    //Actually, I'm reasonably sure this never happens so long as every element is present on every
                    //tree. Assume that for now.
                }
            }

            //need to find the bundles before reattaching anything

            Map<PartitionedTreeNode, List<PartitionedTreeNode>> bundlesForHooks = new HashMap<>();
            Map<PartitionedTreeNode, Double> newHeights = new HashMap<>();
            for(PartitionedTreeNode hook : hooks){
                PartitionedTreeNode hookParent = (PartitionedTreeNode)hook.getParent();
                double extension = hookParent.getHeight() - tree.thisTreeHeight(iP.getHeight());
                double newHeight = tree.thisTreeHeight(newTTheight) + extension;
                newHeights.put(hook, newHeight);

                //at the moment the branch is in the same element as j

                currentBranchElement = j.getPartitionElementNumber();

                while(newHeight > tree.thisTreeHeight(tt.getInfectionHeightByNr(currentBranchElement))){
                    currentBranchElement = tt.getAncestorPartitionElement(currentBranchElement);
                }

                List<PartitionedTreeNode> bundle = new ArrayList<>();

                hrAdjust += Math.log(tree.getBundle((PartitionedTreeNode)tree.getRoot(), extension,
                        currentBranchElement, bundle));

                bundlesForHooks.put(hook, bundle);

            }

            for(PartitionedTreeNode hook : hooks){
                PartitionedTreeNode hookParent = (PartitionedTreeNode)hook.getParent();
                if(hookParent != null){
                    double newHeight = newHeights.get(hook);
                    List<PartitionedTreeNode> bundle = bundlesForHooks.get(hook);

                    PartitionedTreeNode newChild = bundle.get(Randomizer.nextInt(bundle.size()));
                    PartitionedTreeNode newParent = (PartitionedTreeNode)newChild.getParent();

                    if(newParent != null){
                        newParent.removeChild(newChild);
                        newParent.addChild(hookParent);
                        hookParent.addChild(newChild);
                        hookParent.setHeight(newHeight);
                    } else {
                        hookParent.addChild(newChild);
                        tree.setRoot(hookParent);
                    }

                    //The bundle may need to change because of the reinsertion; the branch at that height may have
                    //changed
                    bundlesForHooks.remove(hook);

                    for(PartitionedTreeNode remainingHook : bundlesForHooks.keySet()){
                        List<PartitionedTreeNode> remainingBundle = bundlesForHooks.get(remainingHook);
                        for(PartitionedTreeNode twig : remainingBundle){
                            if(twig == newChild && newHeights.get(remainingHook) > hookParent.getHeight()){
                                remainingBundle.set(remainingBundle.indexOf(twig), hookParent);
                            }
                        }
                    }

                } else {
                    //again, shouldn't happen
                }

            }


        }

        if (jP == null) {
            // adding the node to the root of the tree
            iP.addChild(j);
            tt.setRoot(iP);
        } else {
            // remove destination edge j from its parent
            jP.removeChild(j);

            // add destination edge to the parent of node
            iP.addChild(j);

            // and add the parent of i as a child of the former parent of j
            jP.addChild(iP);
        }

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


        for(GuidedPartitionedTree tree : phylogenies){
            tree.updatePartitions();
        }

        return hrAdjust;

    }

}


