/*
* File GuidedPhylogenySubtreeLeap.java
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
import beast.evolution.datatype.IntegerData;
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
 * This performs the subtree leap move on the single phylogeny from a locus without changing the transmission tree
 */

public class GuidedSubtreeLeapA extends TreeOperator {

    public Input<PartitionedTree> ttInput = new Input<>("tt", "The transmission tree (a partitioned phylogeny with" +
            " second-type rules) to be used as the guide");
    public Input<Double> sizeInput = new Input<>("size", "");

    private PartitionedTree tt;
    private double size;


    public void initAndValidate(){
        tt = ttInput.get();
        if(tt.getRules() != PartitionedTree.Rules.SECOND_TYPE){
            throw new IllegalArgumentException("The guide tree must be partitioned according to second-type rules");
        }
        if(!(treeInput.get() instanceof GuidedPartitionedTree)){
            throw new IllegalArgumentException("This operator is for phylogenies guided by a transmission tree object");
        }

        size = sizeInput.get();
    }

    public double proposal() {

        //this move works only within the transmission chain leading to some clinical case

        int endOfChain = Randomizer.nextInt(tt.getNElements());
        HashMap<Node, Node> references = new HashMap<>();
        GuidedPartitionedTree originalTree = (GuidedPartitionedTree)treeInput.get();
        Tree prunedSubtree = originalTree.extractAncestralSubtree(endOfChain, references);

        //todo are heights actually recalculated?
        double differenceInRootHeights = originalTree.getRoot().getHeight() - prunedSubtree.getRoot().getHeight();
        if(differenceInRootHeights < 0){
            throw new RuntimeException("Pruning must have failed");
        }

        double logq;

        final double delta = getDelta();

        final Node root = prunedSubtree.getRoot();

        Node node;

        do {
            // choose a random node avoiding root. I _think_ that nodes that don't exist in the original tree are
            // OK, since it's their _parents_ that get moved.
            node = prunedSubtree.getNode(Randomizer.nextInt(prunedSubtree.getNodeCount()));

        } while (node == root);

        // get its parent - this is the node we will prune/graft
        final Node parent = node.getParent();

        // get the node's sibling
        final Node sibling = getOtherChild(parent, node);

        // and its grand parent
        final Node grandParent = parent.getParent();

        final Map<Node, Double> destinations = getDestinations(node, parent, sibling, delta);
        final List<Node> destinationNodes = new ArrayList<>(destinations.keySet());

        // pick uniformly from this list
        int r = Randomizer.nextInt(destinations.size());

        double forwardProbability = 1.0 / destinations.size();

        final Node j = destinationNodes.get(r);
        final double newHeight = destinations.get(j);

        final Node jParent = j.getParent();

        if (jParent != null && newHeight > jParent.getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < j.getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        //all operations have to be done to both the copy (for HR purposes) and the original tree. This is annoying.

        Node originalParent = references.get(parent);
        Node originalSibling = references.get(sibling);
        Node originalGrandparent = grandParent != null ? references.get(grandParent) : null;
        Node originalJ = references.get(j);
        Node originalJParent = jParent != null ? references.get(jParent) : null;

        if (j == parent || jParent == parent) {
            // the subtree is not actually moving but the height will change
        } else {
            if (grandParent == null) {
                // if the parent of the original node is the root then the sibling becomes
                // the root.
                parent.removeChild(sibling);
                sibling.setParent(null);
                prunedSubtree.setRoot(sibling);

                originalParent.removeChild(originalSibling);
                originalSibling.setParent(null);
                originalTree.setRoot(originalSibling);

            } else {
                // remove the parent of node by connecting its sibling to its grandparent.
                parent.removeChild(sibling);
                grandParent.removeChild(parent);
                grandParent.addChild(sibling);

                originalParent.removeChild(originalSibling);
                originalGrandparent.removeChild(originalParent);
                originalGrandparent.addChild(originalSibling);
            }

            if (jParent == null) {
                // adding the node to the root of the tree
                parent.addChild(j);
                prunedSubtree.setRoot(parent);

                originalParent.addChild(originalJ);
                originalTree.setRoot(originalParent);
            } else {
                // remove destination edge j from its parent
                jParent.removeChild(j);

                // add destination edge to the parent of node
                parent.addChild(j);

                // and add the parent of i as a child of the former parent of j
                jParent.addChild(parent);

                originalJParent.removeChild(originalJ);
                originalParent.addChild(originalJ);
                originalJParent.addChild(originalParent);
            }
        }

        parent.setHeight(newHeight);
        originalParent.setHeight(newHeight + differenceInRootHeights);

        if (parent.getParent() != null && newHeight > parent.getParent().getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < node.getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < getOtherChild(parent, node).getHeight()) {
            throw new IllegalArgumentException("height error");
        }

        final Map<Node, Double> reverseDestinations = getDestinations(node, parent, getOtherChild(parent, node), delta);
        double reverseProbability = 1.0 / reverseDestinations.size();

        // hastings ratio = reverse Prob / forward Prob
        logq = Math.log(reverseProbability) - Math.log(forwardProbability);
        return logq;
    }

    private Map<Node, Double> getDestinations(Node node, Node parent, Node sibling, double delta) {

        final Map<Node, Double> destinations = new HashMap<>();

        // get the parent's height
        final double height = parent.getHeight();

        final double heightBelow = height - delta;

        if (heightBelow > node.getHeight()) {
            // the destination height below the parent is compatible with the node
            // see if there are any destinations on the sibling's branch
            final List<Node> edges = new ArrayList<Node>();

            getIntersectingEdges(sibling, heightBelow, edges);

            // add the intersecting edges and the height
            for (Node n : edges) {
                destinations.put(n, heightBelow);
            }
        }

        final double heightAbove = height + delta;

        Node node1 = parent;

        // walk up to root
        boolean done = false;
        while (!done) {
            Node parent1 = node1.getParent();

            if (parent1 != null) {
                final double height1 = parent1.getHeight();
                if (height1 < heightAbove) {
                    // haven't reached the height above the original height so go down
                    // the sibling subtree
                    Node sibling1 = getOtherChild(parent1, node1);

                    double heightBelow1 = height1 - (heightAbove - height1);

                    if (heightBelow1 > node.getHeight()) {

                        final List<Node> edges = new ArrayList<Node>();

                        getIntersectingEdges(sibling1, heightBelow1, edges);

                        // add the intersecting edges and the height
                        for (Node n : edges) {
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

    private int getIntersectingEdges(Node node, double height, List<Node> edges) {

        final Node parent = node.getParent();

        if (parent.getHeight() < height) return 0;

        if (node.getHeight() < height) {
            edges.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i <  node.getChildCount(); i++) {
            count += getIntersectingEdges(node.getChild(i), height, edges);
        }
        return count;
    }


}
