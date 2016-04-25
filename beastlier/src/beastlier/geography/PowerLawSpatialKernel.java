/*
* File PowerLawSpatialKernel.java
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

import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class PowerLawSpatialKernel extends SpatialKernel {

    public Input<RealParameter> alphaInput = new Input<>("alpha", "the alpha parameter of this kernel");

    RealParameter alpha;

    public void initAndValidate() {
        alpha = alphaInput.get();
    }

    @Override
    public double getValue(double distance) {
        double alphaValue = alpha.getValue();

        return Math.pow(distance, -alphaValue);
    }
}
