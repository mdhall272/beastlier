/*
* File NormalGamma.java
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
package beast.math.distributions;

import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.GammaFunction;
import org.apache.commons.math.distribution.Distribution;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class NormalGamma extends ParametricDistribution {

    final public Input<RealParameter> muInput = new Input<>("mu", "mu parameter of distribution");
    final public Input<RealParameter> lambdaInput = new Input<>("lambda", "mu parameter of distribution");
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "mu parameter of distribution");
    final public Input<RealParameter> betaInput = new Input<>("beta", "mu parameter of distribution");

    public void initAndValidate() {
    }

    @Override
    public Distribution getDistribution() {
        return null;
    }

    @Override
    public double calcLogP(Function fun) {
        if(fun.getDimension()!=2){
            throw new IllegalArgumentException("Argument to calcLogP has the wrong dimension (should be 2)");
        }

        double mu = getMu();
        double lambda = getLambda();
        double alpha = getAlpha();
        double beta = getBeta();

        double x = fun.getArrayValue(0);
        double tau = fun.getArrayValue(1);

        return alpha*Math.log(beta) + 0.5*Math.log(lambda) + (alpha-0.5)*Math.log(tau)-beta*tau
                -lambda*tau*Math.pow(x-mu,2)/2 - GammaFunction.lnGamma(alpha) - 0.5*Math.log(2*Math.PI);
    }

    public double getMu(){
        return muInput.get().getValue();
    }

    public double getLambda(){
        return lambdaInput.get().getValue();
    }

    public double getAlpha(){
        return alphaInput.get().getValue();
    }

    public double getBeta(){
        return betaInput.get().getValue();
    }

}
