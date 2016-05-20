/*
* File GuidedPartitionedTree.java
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
package beast.evolution.tree;

import beast.core.Input;
import beastlier.outbreak.ClinicalCase;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import java.util.*;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 *
 * A partitioned tree class intended for individual locus phylogenies in a multilocus model. The phylogeny conforms
 * to a transmission structure given by a separate tree with second-type rules, which is the transmission tree.
 */

public class GuidedPartitionedTree extends PartitionedTree {

    public Input<EpidemiologicalPartitionedTree> ttInput = new Input<>("tt", "The transmission tree");

    private EpidemiologicalPartitionedTree tt;

    //what must be added to the height in this tree to get the height in the guide tree
    private double heightAdjustment;

    public void initAndValidate(){
        super.initAndValidate();

        tt = ttInput.get();
        if(tt.getRules() != Rules.SECOND_TYPE){
            throw new IllegalArgumentException("The guide tree must have second-type rules");
        }

        double latestGuideDate = Double.NEGATIVE_INFINITY;

        for(Node node : tt.getExternalNodes()){
            if(tt.getDate(node.getHeight())>latestGuideDate){
                latestGuideDate = tt.getDate(node.getHeight());
            }
        }

        double latestDateHere = Double.NEGATIVE_INFINITY;

        for(Node node : getExternalNodes()){
            if(getDate(node.height)>latestDateHere){
                latestDateHere = getDate(node.height);
            }
        }

        //This is a bad thing that shouldn't happen
        if(latestDateHere > latestGuideDate){
            throw new RuntimeException("Tips in the guide tree cannot be dated earlier than the tips from this tree");
        }

        heightAdjustment = latestGuideDate - latestDateHere;
    }

    public double guideTreeHeight(double thisTreeHeight){
        return thisTreeHeight + heightAdjustment;
    }

    public double thisTreeHeight(double guideTreeHeight){
        return guideTreeHeight - heightAdjustment;
    }

    //Returns the element that the ancestor of this node was present in at the given height (in either this tree or the
    //guide tree.

    private int elementAtHeight(PartitionedTreeNode node, double height, boolean inGuideTree){

        double gtHeight = inGuideTree ? height : guideTreeHeight(height);
        double gtNodeHeight = guideTreeHeight(node.getHeight());

        if(gtHeight < gtNodeHeight){
            throw new IllegalArgumentException("This function is ambiguous if the given height is lower than the " +
                    "node height");
        }

        int nodeElement = node.getPartitionElementNumber();
        PartitionedTreeNode tipInGuideTree = tt.getTip(nodeElement);

        //in the guide tree, internal nodes are transmissions

        PartitionedTreeNode currentNode = tipInGuideTree;
        while(currentNode != null & currentNode.getHeight() < gtHeight){
            currentNode = (PartitionedTreeNode) currentNode.getParent();
        }

        //If the node immediately after this time point is A infecting B, we want to return A
        if(currentNode!=null){
            return currentNode.getPartitionElementNumber();
        }

        //if this height is earlier than the root of the guide tree, need to check the length of its root branch
        //todo At some point this should be a bone fide node, but not today

        if(gtHeight > tt.getRoot().getHeight()+tt.getRootBranchLength() ){
            return -1;
        } else {
            return ((PartitionedTreeNode)tt.getRoot()).getPartitionElementNumber();
        }

    }

    //Partition the internal nodes according to the guide; return false if you can't do it

    public boolean updatePartitions(){

        //this is going to be slow if done like this, but it will catch problems.

        for(Node node : getInternalNodes()){
            PartitionedTreeNode castNode = (PartitionedTreeNode)node;
            List<Integer> requiredElements = new ArrayList<>();

            for(Node tip : node.getAllLeafNodes()){
                PartitionedTreeNode castTip = (PartitionedTreeNode)tip;
                int elementAtNodeHeight = elementAtHeight(castTip, node.getHeight(), false);
                requiredElements.add(elementAtNodeHeight);
            }
            if(requiredElements.size()!=1){
                return false;
            }
            int elementNo = requiredElements.get(0);

            castNode.setPartitionElementNumber(elementNo);

            if(elementNo!=-1) {
                List<Integer> ancestralChain = tt.getAncestralChain(elementNo);

                List<Double> transitionHeights = new ArrayList<>();
                List<Integer> infectorsAlongBranch = new ArrayList<>();

                //starting at the end of the branch


                if (!node.isRoot()) {
                    int currentPositionInChain = 0;
                    double currentHeight = node.getHeight();
                    double parentHeight = node.getParent().getHeight();
                    while (currentHeight < parentHeight) {
                        // going up the branch from the bottom, the arrays list the heights of the next transition (or
                        // the parent height) and the element number in the interval ending in that height

                        currentHeight = tt.getInfectionHeightByNr(ancestralChain.get(currentPositionInChain));
                        if (currentHeight < parentHeight) {
                            transitionHeights.add(currentHeight);
                            infectorsAlongBranch.add(ancestralChain.get(currentPositionInChain));
                            currentPositionInChain++;
                        } else {
                            transitionHeights.add(parentHeight);
                            infectorsAlongBranch.add(ancestralChain.get(currentPositionInChain));
                        }
                    }
                } else {
                    //the root branch stretches to infinity
                    for (int no : ancestralChain) {
                        infectorsAlongBranch.add(no);
                        transitionHeights.add(tt.getInfectionHeightByNr(no));
                    }
                }
                int[] branchElements = Ints.toArray(infectorsAlongBranch);
                double[] branchTransitionHeights = Doubles.toArray(transitionHeights);

                node.setMetaData("historyElements", branchElements);
                node.setMetaData("historyHeights", branchTransitionHeights);
            } else {
                //we're above the root of the guide tree
                List<Double> transitionHeights = new ArrayList<>();
                List<Integer> infectorsAlongBranch = new ArrayList<>();
                transitionHeights.add(node.isRoot() ? Double.POSITIVE_INFINITY : node.getParent().getHeight());
                infectorsAlongBranch.add(-1);

                int[] branchElements = Ints.toArray(infectorsAlongBranch);
                double[] branchTransitionHeights = Doubles.toArray(transitionHeights);

                node.setMetaData("historyElements", branchElements);
                node.setMetaData("historyHeights", branchTransitionHeights);
            }

        }
        return true;
    }


    public boolean somethingIsDirty(){
        return super.somethingIsDirty() || tt.somethingIsDirty();
    }

    //extracts the subtree taking place in and element and all its ancestors in the transmission tree. Probably
    //extremely slow.

    public PartitionedTree extractAncestralSubtree(int elementNo, HashMap<Node, Node> references){

        List<Integer> ancestralChain = getAncestralChain(elementNo);

        PartitionedTreeNode root = (PartitionedTreeNode)getRoot();

        return new PartitionedTree(copyDown(root, references, ancestralChain));
    }

    private PartitionedTreeNode copyDown(PartitionedTreeNode oldNode, HashMap<Node, Node> references,
                                         List<Integer> chain){
        if(chain.contains(oldNode.getPartitionElementNumber())){
            PartitionedTreeNode copyNode = new PartitionedTreeNode();
            copyNode.height = oldNode.getHeight();
            copyNode.labelNr = oldNode.getNr();
            copyNode.metaDataString = oldNode.metaDataString;
            copyNode.setParent(null);
            copyNode.setID(oldNode.getID());
            copyNode.setPartitionElementNumber(oldNode.getPartitionElementNumber());
            references.put(copyNode, oldNode);

            if(!oldNode.isLeaf()){
                for(int i=0; i<oldNode.getChildCount(); i++){
                    copyNode.addChild(copyDown((PartitionedTreeNode)oldNode.getChild(i), references, chain));
                }
            }

            return copyNode;
        } else {
            //put in a new tip to represent the point at which an infection went outside the chain
            int infectedElement = oldNode.getPartitionElementNumber();
            double newHeight = thisTreeHeight(tt.getInfectionHeightByNr(infectedElement));

            PartitionedTreeNode transNode = new PartitionedTreeNode();
            transNode.height = newHeight;
            transNode.setParent(null);
            transNode.setPartitionElementNumber(infectedElement);
            return transNode;
        }
    }


    public void explodeTree(){
        elementsAsTrees = new HashMap<>();
        elementsAsTrees.put(-1, new ArrayList<>());
        for(int i=0; i<getNElements(); i++){
            elementsAsTrees.put(i, new ArrayList<>());
        }
        scanForTreelets((PartitionedTreeNode)getRoot(), elementsAsTrees);
    }

    private Node scanForTreelets(PartitionedTreeNode node, AbstractMap<Integer, List<Treelet>> elementsAsTrees) {

        double[] changeHeights = (double[]) node.getMetaData("historyHeights");
        int[] changeElements = (int[]) node.getMetaData("historyElements");

        if(changeHeights.length == 1){
            //no infections along the branch
            Node copy;
            if(node.isLeaf()){
                copy = new Node(node.getID());
            } else {
                copy = new Node();
                for(int childNo = 0; childNo < node.getChildCount(); childNo++){
                    node.addChild(scanForTreelets((PartitionedTreeNode)node.getChild(childNo), elementsAsTrees));
                }
            }
            copy.setHeight(node.getHeight());
            return copy;
        } else {
            //the section nearest the child
            Node copyRoot;
            if(node.isLeaf()){
                copyRoot = new Node(node.getID());
            } else {
                copyRoot = new Node();
                for(int childNo = 0; childNo < node.getChildCount(); childNo++){
                    node.addChild(scanForTreelets((PartitionedTreeNode)node.getChild(childNo), elementsAsTrees));
                }
            }
            copyRoot.setHeight(node.getHeight());
            Treelet treelet = new Treelet(new Tree(copyRoot), changeHeights[0]-node.getHeight());
            List<Treelet> elementTreelets = elementsAsTrees.get(node.getPartitionElementNumber());
            elementTreelets.add(treelet);

            if(changeElements.length > 2){
                //one or more little segments
                for(int segmentNo = 1; segmentNo < changeElements.length-2; segmentNo++){
                    Node segmentEnd = new Node("Transmission_" + getElementString(changeElements[segmentNo-1]));
                    Treelet treelet1 = new Treelet(new Tree(segmentEnd),
                            changeHeights[segmentNo] - changeHeights[segmentNo-1]);
                    elementTreelets = elementsAsTrees.get(changeElements[segmentNo]);
                    elementTreelets.add(treelet1);
                }
            }

            //the section nearest the parent
            Node transmissionTip = new Node("Transmission_"
                    + getElementString(changeElements[changeElements.length-2]));
            transmissionTip.setHeight(changeHeights[changeElements.length]-2);
            return transmissionTip;
        }

    }


    // The bundle is the set of branches that are in the given element at the given height (even if neither node is
    // in it); more efficient to borrow the intersectingEdges trick from operators, which also helps if some subtrees
    // are temporarily disconnected mid-proposal

    public int getBundle(PartitionedTreeNode node, double height, int elementNo, List<PartitionedTreeNode> twigs){

        List<Integer> ancestralChain = getAncestralChain(elementNo);

        final PartitionedTreeNode parent = (PartitionedTreeNode)node.getParent();

        //we can stop if:
        //1) We're so far down that the edge cannot intersect
        //2) We're higher up than the required height and in a partition element which is not part of the ancestral
        //chain. There's no way a lineage can re-enter the required element subsequently.

        if((parent != null && parent.getHeight() < height) ||
                (node.getHeight() > height && !(ancestralChain.contains(node.getPartitionElementNumber())))){
            return 0;
        }

        if(parent!=null && !(ancestralChain.contains(parent.getPartitionElementNumber()))){
            throw new RuntimeException("The bundle function is not working");
        }

        if (node.getHeight() < height) {
            if(parent == null || ancestralChain.contains(parent.getPartitionElementNumber())){
                int[] branchElements = (int[])node.getMetaData("historyElements");
                double[] transitionHeights = (double[])node.getMetaData("historyHeights");

                double currentHeight = node.getHeight();
                int currentIndex = -1;
                int currentElement = node.getPartitionElementNumber();

                while(currentHeight < height){
                    currentIndex ++;
                    currentHeight = transitionHeights[currentIndex];
                    currentElement = branchElements[currentIndex];
                }

                if(currentElement == elementNo){
                    if(twigs != null) {
                        twigs.add(node);
                    }
                    return 1;
                }

            }
        }

        int count = 0;
        for (int i = 0; i < node.getChildCount(); i++) {
            count += getBundle((PartitionedTreeNode)node.getChild(i), height, elementNo, twigs);
        }
        return count;
    }

    // The hooks are where the phylogeny nodes protruding from the transmission tree branch as it is moved hook
    // onto the rest of the tree. This can feasibly be the root's "parent", and because of difficulties with
    // the root otherwise, this function returns the ends of the hook branches, not the beginnings.
    // They can be found by finding the bundle at the point of infection and then tracing up to the nodes below
    // the points of attachment to the rest of the tree

    public List<PartitionedTreeNode> getHooks(int elementNo, double height){
        // First, get all the nodes whose parent branches are at in the required element at the right time (even if
        // neither node is

        List<PartitionedTreeNode> laterNodes = new ArrayList<>();

        getBundle((PartitionedTreeNode)getRoot(), height, elementNo, laterNodes);

        //but this isn't quite it; if any two nodes in that set coalesce with each other before they coalesce
        //which nodes from elsewhere in the tree, then the hook extends from their ancestor, not them.

        //dear God is this method inefficient right now

        boolean satisfied = false;

        while(!satisfied){
            Set<PartitionedTreeNode> parents = new HashSet<>();
            boolean rootIsIn = false;
            for(PartitionedTreeNode child : laterNodes){
                if(child.getParent() != null){
                    parents.add((PartitionedTreeNode)child.getParent());
                } else {
                    rootIsIn = true;
                }
            }
            //if all the parents are distinct then you can stop.
            if((rootIsIn && parents.size() == laterNodes.size()-1)
                    || (!(rootIsIn) && parents.size() == laterNodes.size())){
                satisfied = true;
            } else {
                outer_loop:
                for(PartitionedTreeNode node1 : laterNodes){
                    for(PartitionedTreeNode node2 : laterNodes){
                        if(node1!=node2){
                            if(node1.getParent() == node2.getParent()){
                                laterNodes.remove(node1);
                                laterNodes.remove(node2);
                                laterNodes.add((PartitionedTreeNode)node1.getParent());
                                break outer_loop;
                            }
                        }
                    }
                }
            }
        }
        return laterNodes;
    }
}




