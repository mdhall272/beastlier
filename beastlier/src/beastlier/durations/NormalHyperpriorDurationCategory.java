/*
* File NormalHyperpriorDurationCategory.java
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
import beast.math.GammaFunction;
import beast.math.distributions.NormalGamma;
import beastlier.outbreak.*;

import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("A group of clinical cases whose latent or infectious periods are assumed to be drawn from the same" +
        " unknown normal distribution")
public class NormalHyperpriorDurationCategory extends DurationCategory {

    public Input<NormalGamma> distributionInput = new Input<>("NormalGamma", "The normal-gamma distribution " +
            "determining the hyperprior");

    private NormalGamma distribution;

    @Override
    public void initAndValidate() {
        distribution = distributionInput.get();
        hasProbability = true;
    }

    //todo posterior distribution of these parameters

    public double getLogProbability(Function values){

        int count = values.getDimension();

        double mu_0 = distribution.getMu();
        double lambda_0 = distribution.getLambda();
        double alpha_0 = distribution.getAlpha();
        double beta_0 = distribution.getBeta();

        double lambda_n = lambda_0 + count;
        double alpha_n = alpha_0 + count/2;
        double sum = 0;
        for (Double infPeriod : values.getDoubleValues()) {
            sum += infPeriod;
        }
        double mean = sum/count;

        double sumOfDifferences = 0;
        for (Double duration : values.getDoubleValues()) {
            sumOfDifferences += Math.pow(duration-mean,2);
        }


        double beta_n = beta_0 + 0.5*sumOfDifferences
                + lambda_0*count*Math.pow(mean-mu_0, 2)/(2*(lambda_0+count));


        return GammaFunction.lnGamma(alpha_n)
                - GammaFunction.lnGamma(alpha_0)
                + alpha_0*Math.log(beta_0)
                - alpha_n*Math.log(beta_n)
                + 0.5*Math.log(lambda_0)
                - 0.5*Math.log(lambda_n)
                - (count/2)*Math.log(2*Math.PI);

    }
}
