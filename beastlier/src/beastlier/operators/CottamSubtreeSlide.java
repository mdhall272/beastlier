/*
* File SecondTypeSubtreeSlide.java
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

import beast.evolution.operators.SubtreeSlide;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 *
 * The base operator was clearly not intended to be subclassed...
 */

public class CottamSubtreeSlide extends SubtreeSlide {

    double size;
    double limit;

    public void initAndValidate() {
        size = sizeInput.get();
        limit = limitInput.get();
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules == PartitionedTree.Rules.COTTAM)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the second type");
        }
    }

    @Override
    public double proposal() {
        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        double logq;

        PartitionedTreeNode i;
        final boolean markClades = markCladesInput.get();
        // 1. choose a random node avoiding root
        final int nodeCount = tree.getNodeCount();
        do {
            i = (PartitionedTreeNode)tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot());

        final PartitionedTreeNode p = (PartitionedTreeNode)i.getParent();
        final PartitionedTreeNode CiP = (PartitionedTreeNode)getOtherChild(p, i);
        final PartitionedTreeNode PiP = (PartitionedTreeNode)p.getParent();

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = p.getHeight();
        final double newHeight = oldHeight + delta;

        // 3. if the move is up
        if (delta > 0) {

            // 3.1 if the topology will change
            if (PiP != null && PiP.getHeight() < newHeight) {
                // find new parent
                PartitionedTreeNode newParent = PiP;
                PartitionedTreeNode newChild = p;
                while (newParent.getHeight() < newHeight) {
                    newChild = newParent;
                    if( markClades ) newParent.makeDirty(Tree.IS_FILTHY); // JH
                    newParent = (PartitionedTreeNode)newParent.getParent();
                    if (newParent == null) break;
                }
                // the moved node 'p' would become a child of 'newParent'
                //

                // 3.1.1 if creating a new root
                if (newChild.isRoot()) {
                    replace(PiP, p, CiP);
                    //p's subtree is now disconnected from the tree
                    if(i.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                        PartitionedTreeNode currentNode = PiP;
                        while (currentNode!=null &&
                                currentNode.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                            currentNode.setPartitionElementNumber(CiP.getPartitionElementNumber());
                            currentNode.setPartitionDirty(true);
                            currentNode = (PartitionedTreeNode) currentNode.getParent();
                        }
                    }

                    replace(p, CiP, newChild);
                    p.setParent(null);
                    tree.setRoot(p);
                    // and it's back. No need to worry about anything further up the tree since it's now the root.
                    if(i.getPartitionElementNumber() != p.getPartitionElementNumber()){
                        p.setPartitionElementNumber(newChild.getPartitionElementNumber());
                        p.setPartitionDirty(true);
                    }
                }
                // 3.1.2 no new root
                else {
                    replace(PiP, p, CiP);
                    //gone
                    if(i.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                        PartitionedTreeNode currentNode = PiP;
                        while (currentNode != null &&
                                currentNode.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                            currentNode.setPartitionElementNumber(CiP.getPartitionElementNumber());
                            currentNode.setPartitionDirty(true);
                            currentNode = (PartitionedTreeNode) currentNode.getParent();
                        }
                    }

                    replace(p, CiP, newChild);
                    replace(newParent, newChild, p);
                    //back
                    if(i.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                        PartitionedTreeNode currentNode = newParent;
                        while (currentNode != null &&
                                currentNode.getPartitionElementNumber() == newChild.getPartitionElementNumber()) {
                            currentNode.setPartitionElementNumber(p.getPartitionElementNumber());
                            currentNode.setPartitionDirty(true);
                            currentNode = (PartitionedTreeNode) currentNode.getParent();
                        }
                    } else {
                        p.setPartitionElementNumber(newChild.getPartitionElementNumber());
                        p.setPartitionDirty(true);
                    }
                }

                p.setHeight(newHeight);

                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChild, oldHeight, null);
                //System.out.println("possible sources = " + possibleSources);

                logq = -Math.log(possibleSources);

            } else {
                // just change the node height
                p.setHeight(newHeight);
                logq = 0.0;
            }
        }
        // 4 if we are sliding the subtree down.
        else {

            // 4.0 is it a valid move?
            if (i.getHeight() > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (CiP.getHeight() > newHeight) {

                final List<Node> newChildren = new ArrayList<>();
                final int possibleDestinations = intersectingEdges(CiP, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = Randomizer.nextInt(newChildren.size());
                final PartitionedTreeNode newChild = (PartitionedTreeNode)newChildren.get(childIndex);
                final PartitionedTreeNode newParent = (PartitionedTreeNode)newChild.getParent();

                // 4.1.1 if p was root
                if (p.isRoot()) {
                    // new root is CiP
                    // nothing to do when removing it as root

                    replace(p, CiP, newChild);
                    replace(newParent, newChild, p);

                    CiP.setParent(null);
                    tree.setRoot(CiP);

                    if(i.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                        PartitionedTreeNode currentNode = newParent;
                        while (currentNode != null &&
                                currentNode.getPartitionElementNumber() == newChild.getPartitionElementNumber()) {
                            currentNode.setPartitionElementNumber(p.getPartitionElementNumber());
                            currentNode.setPartitionDirty(true);
                            currentNode = (PartitionedTreeNode) currentNode.getParent();
                        }
                    } else {
                        p.setPartitionElementNumber(newChild.getPartitionElementNumber());
                        p.setPartitionDirty(true);
                    }

                } else {
                    replace(PiP, p, CiP);
                    //gone
                    if(i.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                        PartitionedTreeNode currentNode = PiP;
                        while (currentNode != null &&
                                currentNode.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                            currentNode.setPartitionElementNumber(CiP.getPartitionElementNumber());
                            currentNode.setPartitionDirty(true);
                            currentNode = (PartitionedTreeNode) currentNode.getParent();
                        }
                    }

                    replace(p, CiP, newChild);
                    replace(newParent, newChild, p);
                    //back
                    if(i.getPartitionElementNumber() == p.getPartitionElementNumber()) {
                        PartitionedTreeNode currentNode = newParent;
                        while (currentNode != null &&
                                currentNode.getPartitionElementNumber() == newChild.getPartitionElementNumber()) {
                            currentNode.setPartitionElementNumber(p.getPartitionElementNumber());
                            currentNode.setPartitionDirty(true);
                            currentNode = (PartitionedTreeNode) currentNode.getParent();
                        }
                    } else {
                        p.setPartitionElementNumber(newChild.getPartitionElementNumber());
                        p.setPartitionDirty(true);
                    }
                }

                p.setHeight(newHeight);
                if( markClades ) {
                    // make dirty the path from the (down) moved node back up to former parent.
                    Node n = p;
                    while( n != CiP ) {
                        n.makeDirty(Tree.IS_FILTHY); // JH
                        n = n.getParent();
                    }
                }

                logq = Math.log(possibleDestinations);
            } else {
                p.setHeight(newHeight);
                logq = 0.0;
            }
        }

        return logq;
    }

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

}
