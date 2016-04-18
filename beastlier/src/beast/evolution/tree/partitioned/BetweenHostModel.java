/*
* File BetweenHostModel.java
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
package beast.evolution.tree.partitioned;

import beast.core.Input;
import beast.evolution.tree.EpidemiologicalPartitionedTree;
import beast.evolution.tree.PartitionedTree;
import beast.evolution.tree.TreeDistribution;
import beastlier.outbreak.Outbreak;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public abstract class BetweenHostModel extends TreeDistribution {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The collection of clinical cases");

    private EpidemiologicalPartitionedTree tree;
    private Outbreak outbreak;

    public void initAndValidate(){
        if(!(treeInput.get() instanceof EpidemiologicalPartitionedTree)){
            throw new IllegalArgumentException("Trees given to the between-host model must have node partitions and" +
                    " an outbreak");
        }

        tree = (EpidemiologicalPartitionedTree) treeInput.get();
        outbreak = outbreakInput.get();
    }



}
