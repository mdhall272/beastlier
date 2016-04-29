/*
* File TransmissionTreeStuctureLogger.java
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
import beast.evolution.tree.PartitionedTree;

import java.io.PrintStream;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class TransmissionTreeStructureLogger extends BEASTObject implements Loggable {

    public Input<PartitionedTree> partitionedTreeInput = new Input<>("partitionedTree", "Node-partitioned tree to log.",
            Input.Validate.REQUIRED);

    PartitionedTree pTree;

    @Override
    public void initAndValidate() {
        pTree = partitionedTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        for(String elementName : pTree.getElementList()){
            out.print(elementName + "\t" );
        }
    }

    @Override
    public void log(int sample, PrintStream out) {
        for(String elementName : pTree.getElementList()){
            //element number
            String infector = pTree.getAncestorPartitionElement(elementName);
            out.print((infector==null ? "start" : infector) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        //nothing to do
    }
}
