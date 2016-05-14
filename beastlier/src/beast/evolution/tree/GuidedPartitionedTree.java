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

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 *
 * A partitioned tree class intended for individual locus phylogenies in a multilocus model. The phylogeny conforms
 * to a transmission structure given by a separate tree with second-type rules, which is the transmission tree.
 */

public class GuidedPartitionedTree extends PartitionedTree {

    public Input<EpidemiologicalPartitionedTree> ttInput = new Input<>("tt", "The transmission tree");

    private EpidemiologicalPartitionedTree tt;


    //Returns the element that the ancestor of this tip was present in at the given height (in either this tree or the
    //guide tree.

    private int elementAtHeight(PartitionedTreeNode node, double height, boolean inGuideTree){
        return 0;
    }

    //Partition the internal nodes according to the guide; return false if you can't do it

    private boolean updatePartitions(){

        return false;

    }

}
