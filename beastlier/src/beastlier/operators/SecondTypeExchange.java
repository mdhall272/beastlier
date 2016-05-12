/*
* File SecondTypeExchange.java
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

import beast.evolution.operators.Exchange;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class SecondTypeExchange extends Exchange {

    @Override
    public void initAndValidate() {
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules == PartitionedTree.Rules.SECOND_TYPE)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the second type");
        }
    }

    protected void exchangeNodes(Node i, Node j,
                                 Node p, Node jP) {
        replace(p, i, j);
        replace(jP, j, i);

        PartitionedTreeNode castI = (PartitionedTreeNode)i;
        PartitionedTreeNode castJ = (PartitionedTreeNode)j;
        PartitionedTreeNode castIP = (PartitionedTreeNode)p;
        PartitionedTreeNode castJP = (PartitionedTreeNode)jP;


        PartitionedTreeNode currentParent = castJP;
        while(currentParent!=null && currentParent.getPartitionElementNumber()==castJ.getPartitionElementNumber()){
            currentParent.setPartitionElementNumber(castI.getPartitionElementNumber());
            currentParent.setPartitionDirty(true);
            currentParent = (PartitionedTreeNode)currentParent.getParent();
        }

        currentParent = castIP;
        while(currentParent!=null && currentParent.getPartitionElementNumber()==castI.getPartitionElementNumber()){
            currentParent.setPartitionElementNumber(castJ.getPartitionElementNumber());
            currentParent.setPartitionDirty(true);
            currentParent = (PartitionedTreeNode)currentParent.getParent();
        }
    }
}

