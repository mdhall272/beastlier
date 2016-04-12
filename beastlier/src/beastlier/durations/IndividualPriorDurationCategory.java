/*
* File IndividualPriorDurationCategory.java
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

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.math.distributions.ParametricDistribution;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("A group of clinical cases whose latent or infectious periods each have the same prior distribution")
public class IndividualPriorDurationCategory extends DurationCategory {

    public Input<ParametricDistribution> distributionInput = new Input<>("distribution", "The prior distribution for" +
            " the length of this duration (latent or infectious period) for all clinical cases in this category");

    private ParametricDistribution distribution;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        distribution = distributionInput.get();
    }



    public double getLogProbability(Function values){
        if(values.getDimension() != cases.size()){
            throw new RuntimeException("Querying the probability of a set of durations of the wrong dimension");
        }

        return distribution.calcLogP(values);
    }

}
