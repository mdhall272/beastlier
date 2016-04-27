/*
* File TimedTransmissionTreeLogger.java
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
package beastlier.util;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.*;

import java.io.PrintStream;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class PartitionedTreeLogger extends BEASTObject implements Loggable {

    public Input<PartitionedTree> partitionedTreeInput = new Input<>(
            "partitionedTree",
            "Node-partitioned tree to log.",
            Input.Validate.REQUIRED);
    public Input<Boolean> fancyInput = new Input<>("fancy", "Output trees with infection events as degree two " +
            "internal nodes", false, Input.Validate.OPTIONAL);

    PartitionedTree pTree;
    boolean fancy;

    @Override
    public void initAndValidate() {
        pTree = partitionedTreeInput.get();
        fancy = fancyInput.get();

        if(fancy && !(pTree instanceof EpidemiologicalPartitionedTree)){
            throw new IllegalArgumentException("Fancy logger only works on EpidemiologicalPartitionedTrees");
        }
    }

    @Override
    public void init(PrintStream out) {
        pTree.init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {

        if(fancy){
            EpidemiologicalPartitionedTree epTree = (EpidemiologicalPartitionedTree)pTree;

            Tree fancyTree = epTree.getFlattenedTree();

            for (Node node : fancyTree.getNodesAsArray()) {
                PartitionedTreeNode pNode = (PartitionedTreeNode) node;
                if(pNode.getPartitionElementNumber()==-1){
                    pNode.metaDataString = pTree.getElementLabel()
                            + "=\""
                            + "start"
                            + "\"";

                } else {
                    pNode.metaDataString = pTree.getElementLabel()
                            + "=\""
                            + pTree.getElementString(pNode.getPartitionElementNumber())
                            + "\"";
                }
            }

            out.print("tree STATE_" + nSample + " = ");
            out.print(fancyTree.getRoot().toSortedNewick(new int[1], true));
            out.print(";");


        } else {

            for (Node node : pTree.getNodesAsArray()) {
                PartitionedTreeNode pNode = (PartitionedTreeNode) node;
                pNode.metaDataString = pTree.getElementLabel()
                        + "=\""
                        + pTree.getElementString(pNode.getPartitionElementNumber())
                        + "\"";
            }

            out.print("tree STATE_" + nSample + " = ");
            out.print(pTree.getRoot().toSortedNewick(new int[1], true));
            out.print(";");
        }
    }

    @Override
    public void close(PrintStream out) {
        pTree.close(out);
    }

}
