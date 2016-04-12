/*
* File DurationCategory.java
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
package beastlier.durations;

import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.math.distributions.ParametricDistribution;
import beastlier.outbreak.ClinicalCase;

import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public abstract class DurationCategory extends BEASTObject {

    //todo some class (between-host model?) should check the categories are mutually exclusive

    public Input<List<ClinicalCase>> casesInput = new Input<>("cases", "A list of cases in this category");

    protected List<ClinicalCase> cases;

    @Override
    public void initAndValidate() {
        cases = casesInput.get();
    }

    public List<ClinicalCase> getCases(){
        return cases;
    }

    // No probability calculations in the base version;

    public double getLogProbability(Function values){
        return 1;
    }
}
