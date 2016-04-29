/*
* File SpatialKernel.java
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

package beastlier.geography;

import beast.core.CalculationNode;
import beast.core.Description;
import beastlier.outbreak.GeographicallyLocatedClinicalCase;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;

import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("A spatial transmission kernel function")
public abstract class SpatialKernel extends CalculationNode implements UnivariateRealFunction {

    public abstract List<String> getParameterIds();

    public double getValue(GeographicallyLocatedClinicalCase case1, GeographicallyLocatedClinicalCase case2) throws
            FunctionEvaluationException{
        return value(Util.EuclideanDistance(case1.getCoordinates(), case2.getCoordinates()));

    }

    public static class Util{

        public static double EuclideanDistance(double[] point1, double[] point2){
            return Math.sqrt(Math.pow(point1[0]-point2[0],2) + Math.pow(point1[1]-point2[1], 2));
        }

    }
}
