/*
* File SecondTypeFlipper.java
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

import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.PartitionedTreeNode;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class SecondTypeFlipper extends TreeOperator {

    public void initAndValidate() {
        if(!(treeInput.get() instanceof PartitionedTree)){
            throw new RuntimeException("This operator is designed for partitioned trees only");
        }
        if(!(((PartitionedTree) treeInput.get()).rules == PartitionedTree.Rules.SECOND_TYPE)){
            throw new RuntimeException("This operator is designed for trees with partition rules of the second type");
        }
    }


    public double proposal() {

        //note that in second type trees there's no need to use the partition-dirty flag, as that is used by the
        //within-host model

        final PartitionedTree tree = (PartitionedTree)treeInput.get(this);

        // find an internal node whose partition we are going to flip from one child to the other
        int nodeToAdjust = Randomizer.nextInt(tree.getInternalNodeCount());
        // if the infection event is the seed of the epidemic, we need to try again

        PartitionedTreeNode parentNode = (PartitionedTreeNode) tree.getInternalNodes().get(nodeToAdjust);

        double hr = adjustTree(parentNode);

        return hr;
    }

    private double adjustTree(PartitionedTreeNode parent) {

        int infectorCase = parent.getPartitionElementNumber();

        //sanity checks

        if(parent.getChildren().size()!=2) {
            throw new RuntimeException("Bifurcating trees only");
        }
        PartitionedTreeNode child1 = (PartitionedTreeNode) parent.getChild(0);
        PartitionedTreeNode child2 = (PartitionedTreeNode) parent.getChild(1);
        if(parent.getPartitionElementNumber() != child1.getPartitionElementNumber() &&
                parent.getPartitionElementNumber() != child2.getPartitionElementNumber()){
            throw new RuntimeException("This tree does not obey second-type partition rules");
        }
        if(child1.getPartitionElementNumber() == child2.getPartitionElementNumber()){
            throw new RuntimeException("This tree does not obey second-type partition rules");
        }

        int infectedCase = parent.getPartitionElementNumber() == child1.getPartitionElementNumber() ?
                child2.getPartitionElementNumber() : child1.getPartitionElementNumber();

        PartitionedTreeNode currentNode = parent;
        while(currentNode!=null && currentNode.getPartitionElementNumber()==infectorCase){
            currentNode.setPartitionDirty(true);
            currentNode.setPartitionElementNumber(infectedCase);
            currentNode = (PartitionedTreeNode) currentNode.getParent();
        }

        child1.setPartitionDirty(true);
        child2.setPartitionDirty(true);

        treeInput.get().setSomethingIsDirty(true);

        //this is a really simple and pleasing move

        return 0;
    }


}
