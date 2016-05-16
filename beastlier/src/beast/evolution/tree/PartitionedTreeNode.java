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

package beast.evolution.tree;

import beast.core.Description;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk> based on MultiTypeNode.java from the MultiTypeTree package by Tim Vaughan
 */
@Description("A node in a beast.evolution.tree with node partitions")
public class PartitionedTreeNode extends Node {

    int partitionElementNumber = -1;
/*
    a node is partition-dirty if something has changed about it that affects the subtree that involves it. This can be:

    1) A change to its height, the height of its parent, or the height of its children. (If a branch connecting
    to parent or child is an infection branch, this will change the subtree even if no nodes _in_ the subtree change
    height.
    2) It has been assigned a new parent or child. This will only not affect the subtree if the new parent/child has
    the same height as the old one, which is sufficiently unlikely as to not be worth considering
    3) It has been moved from one partition element to another
    4) A child has been moved from one partition element to another
    4) Its parent has been moved from _its_ partition element to another or vice versa. Because of how the qis work,
    its parent moving to a from one different element to another is not a problem

    */

    private boolean partitionDirty = true;

    /**
     * Obtain partition element number for this node.
     *
     * @return partition number of this node
     */
    public int getPartitionElementNumber() {
        return partitionElementNumber;
    }

    /**
     * Set partition element number.
     *
     * @param partitionElementNumber New node type.
     */
    public void setPartitionElementNumber(int partitionElementNumber) {
        this.partitionElementNumber = partitionElementNumber;
    }

    /**
     * @return shallow copy of node
     */
    public PartitionedTreeNode shallowCopy() {
        PartitionedTreeNode node = new PartitionedTreeNode();
        node.height = height;
        node.setParent(getParent());
        for(Node child : getChildren()) {
            node.addChild(child);
        }

        node.partitionElementNumber = partitionElementNumber;

        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.ID = ID;

        return node;
    }

    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */


    /**
     * methods from Node.java modified to include type structure
     */


    /**
     * @return (deep) copy of node
     */
    @Override
    public PartitionedTreeNode copy() {
        PartitionedTreeNode node = new PartitionedTreeNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.setParent(null);
        node.ID = ID;
        node.partitionElementNumber = partitionElementNumber;

        if (getLeft()!=null) {
            node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                node.setRight(getRight().copy());
                node.getRight().setParent(node);
            }
        }

        return node;
    }



    /**
     * assign values from a beast.evolution.tree in array representation *
     * @param nodes
     * @param node
     */
    @Override
    public void assignFrom(Node[] nodes, Node node) {
        height = node.getHeight();
        labelNr = node.getNr();
        metaDataString = node.metaDataString;
        setParent(null);
        ID = node.getID();

        if(node instanceof PartitionedTreeNode) {

            PartitionedTreeNode mtNode = (PartitionedTreeNode) node;
            partitionElementNumber = mtNode.partitionElementNumber;
        }
//        if (mtNode.likes==null) {likes = null;}else{
//        	likes = Arrays.copyOf(mtNode.likes,mtNode.likes.length);
//        }

        if (node.getLeft()!=null) {
            setLeft(nodes[node.getLeft().getNr()]);
            getLeft().assignFrom(nodes, node.getLeft());
            getLeft().setParent(this);
            if (node.getRight()!=null) {
                setRight(nodes[node.getRight().getNr()]);
                getRight().assignFrom(nodes, node.getRight());
                getRight().setParent(this);
            }
        }
    }

    void setParent(final Node parent, final boolean inOperator) {
        super.setParent(parent, inOperator);
        setPartitionDirty(true);
        if(parent!=null) {
            ((PartitionedTreeNode) parent).setPartitionDirty(true);
        }
    }

    public void setChild(final int childIndex, final Node node) {
        super.setChild(childIndex, node);
        setPartitionDirty(true);
        ((PartitionedTreeNode)node).setPartitionDirty(true);
    }


    public void setLeft(final Node m_left) {
        super.setLeft(m_left);
        setPartitionDirty(true);
        ((PartitionedTreeNode)m_left).setPartitionDirty(true);
    }


    public void setRight(final Node m_right) {
        super.setRight(m_right);
        setPartitionDirty(true);
        ((PartitionedTreeNode)m_right).setPartitionDirty(true);
    }

    public void setHeight(final double height) {
        super.setHeight(height);
        setPartitionDirty(true);
        PartitionedTreeNode parent = (PartitionedTreeNode)getParent();
        if(parent != null){
            parent.setPartitionDirty(true);
        }
        for(Node child : getChildren()){
            PartitionedTreeNode castChild = (PartitionedTreeNode)child;
            castChild.setPartitionDirty(true);
        }
    }

    /**
     * scale height of this node and all its descendants
     *
     * @param scale scale factor
     */
    public void scale(final double scale) {
        if(getParent()!=null) {
            ((PartitionedTreeNode) getParent()).setPartitionDirty(true);
        }
        setPartitionDirty(true);
        super.scale(scale);
    }

    public void setPartitionDirty(boolean value){
        partitionDirty = value;
    }

    public boolean isPartitionDirty() {
        return partitionDirty;
    }

}