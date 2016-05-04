/*
* File GammaDistributionImplAlt.java
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
package beast.util;

import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.special.Gamma;

/**
 * Simple superclass of GammaDistributionImpl that calculates logP without converting to real space and hence will
 * round densities to zero less readily
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class GammaDistributionImplAlt extends GammaDistributionImpl {

    public GammaDistributionImplAlt(double alpha, double beta) {
        super(alpha, beta);
    }

    public GammaDistributionImplAlt(double alpha, double beta, double inverseCumAccuracy) {
        super(alpha, beta, inverseCumAccuracy);
    }

    public double logDensity(double x) {

        double shape = getAlpha();
        double scale = getBeta();

        return ((shape - 1.0) * (Math.log(x) - Math.log(scale)) - x / scale - Gamma.logGamma(shape)
                - Math.log(scale));
    }

}
