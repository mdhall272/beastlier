/*
 * Copyright (C) 2014 Nicola De Maio
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package tree;


import beast.core.Description;
import beast.evolution.tree.Node;

//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.List;

/**
 *
 * @author Nicola De Maio
 */
@Description("A node that belongs to a partition element.")
public class PElementMemberNode extends Node {

    // Type metadata:
    int pElementNumber = -1;
    
    /**
     * Obtain type value at this node.
     *
     * @return type value at this node
     */
    public int getPElementNumber() {
        return pElementNumber;
    }

    /**
     * Sets type of node.
     *
     * @param pElementNumber New element number.
     */
    public void setpElementNumber(int pElementNumber) {
        //startEditing();
        this.pElementNumber = pElementNumber;
    }

    
    /**
     * Remove all type changes from branch above node.
     */
    public void clearChanges() {
        startEditing();
    }

    /**
     * @return shallow copy of node
     */
    public PElementMemberNode shallowCopy() {
        PElementMemberNode node = new PElementMemberNode();
        node.height = height;
        node.setParent(node.getParent());

        for(Node child : node.getAllChildNodes()){
            node.addChild(child);
        }

        node.pElementNumber = pElementNumber;
                
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
    public PElementMemberNode copy() {
        PElementMemberNode node = new PElementMemberNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.setParent(null);
        node.ID = ID;
        node.pElementNumber = pElementNumber;
//        if (likes==null) {node.likes = null;}else{
//        	node.likes = Arrays.copyOf(likes,likes.length);
//        }
        
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
     * assign values from a tree in array representation *
     * @param nodes
     * @param node
     */
    @Override
    public void assignFrom(Node[] nodes, Node node) {
        height = node.getHeight();
        labelNr = node.getNr();
        metaDataString = node.metaDataString;
        node.setParent(null);
        ID = node.getID();
        
        PElementMemberNode mtNode = (PElementMemberNode) node;
        pElementNumber = mtNode.pElementNumber;
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
}
