/*
* File WithinHostModel.java
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
import beast.evolution.tree.*;
import beastlier.outbreak.Outbreak;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public abstract class WithinHostModel extends TreeDistribution {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The collection of clinical cases");

    protected PartitionedTree tree;
    protected Outbreak outbreak;



    @Override
    public void initAndValidate(){
        if(!(treeInput.get() instanceof EpidemiologicalPartitionedTree)
                && !(treeInput.get() instanceof GuidedPartitionedTree)){
            throw new IllegalArgumentException("Trees given to the within-host model must either have node " +
                    "partitions and an outbreak, or be linked to another one that does");
        }

        tree = (PartitionedTree) treeInput.get();
        outbreak = outbreakInput.get();

        if(tree.getRules() == PartitionedTree.Rules.COTTAM){
            throw new IllegalArgumentException("Trees must be partitioned by third-type rules to have a within-" +
                    "host model");
        }
    }

}
