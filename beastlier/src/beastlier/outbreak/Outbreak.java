/*
* File Outbreak.java
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

package beastlier.outbreak;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.StandardData;

import java.util.Arrays;
import java.util.List;

/**
 * Created by twoseventwo on 11/04/2016.
 */
public abstract class Outbreak extends BEASTObject {

    // todo check you actually still need the taxa here; might not

    public Input<StandardData> caseDataType;
    public Input<TaxonSet> taxa;
    public Input<BooleanParameter> hasLatentPeriods;
    public Input<BooleanParameter> hasGeography;
    public Input<List<ClinicalCase>> cases;


}
