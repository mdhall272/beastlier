/*
* File WithinHostCoalescent.java
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
package beast.evolution.tree.partitioned;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.coalescent.*;
import beast.math.Binomial;
import beastlier.outbreak.ClinicalCase;
import beast.util.BigDecimalUtils;

import java.math.BigDecimal;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("A parameteric within-host coalescent function shared by every clinical case")
public class WithinHostCoalescent extends WithinHostModel {

    public Input<PopulationFunction> functionInput = new Input<>("populationFunction", "The within-host coalescent " +
            "process");

    private PopulationFunction popFunction;
    protected static double tolerance = 1E-10;

    public void initAndValidate(){
        super.initAndValidate();

        popFunction = functionInput.get();
    }


    public double calculateLogP() {

        //checkPartitions();

        double logL = 0;

        for (ClinicalCase aCase : outbreak.getEverInfectedCases()) {

            int number = outbreak.getCaseIndex(aCase);


            //recalculate everything for now

//                if (recalculateCoalescentFlags[number]) {

            Treelet treelet = elementsAsTrees.get(aCase);

            if (treelet.getLeafNodeCount() > 1) {
                TreeIntervals intervals = new TreeIntervals(treelet);

//                    partitionTreeLogLikelihoods[number] = coalescent.calculateLogLikelihood();


                logL += calculateTreeletLogLikelihood(intervals, popFunction, 0, treelet.getZeroHeight());

            } else {
                logL += 0.0;
            }
//            recalculateCoalescentFlags[number] = false;
//        } else {
//            coalescencesLogLikelihood += partitionTreeLogLikelihoods[number];
//        }

//        else {
//                recalculateCoalescentFlags[number] = false;
//            }


        }

        return logL;
    }


    public static double calculateTreeletLogLikelihood(IntervalList intervals, PopulationFunction demographicFunction,
                                                       double threshold, double zeroHeight) {

        double logL = 0.0;

        double startTime = -zeroHeight;
        final int n = intervals.getIntervalCount();

        //TreeIntervals sets up a first zero-length interval with a lineage count of zero - skip this one

        for (int i = 0; i < n; i++) {

            // time zero corresponds to the date of first infection

            final double duration = intervals.getInterval(i);
            final double finishTime = startTime + duration;

            // if this has happened the run is probably pretty unhappy

            if (finishTime == 0) {
                return Double.NEGATIVE_INFINITY;
            }

            final double intervalArea = demographicFunction.getIntegral(startTime, finishTime);
            final double normalisationArea = demographicFunction.getIntegral(startTime, 0);

            if (intervalArea == 0 && duration > tolerance) {
                return Double.NEGATIVE_INFINITY;
            }

            final int lineageCount = intervals.getLineageCount(i);

            if (lineageCount >= 2) {

                final double kChoose2 = Binomial.choose2(lineageCount);

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                    logL += -kChoose2 * intervalArea;

                    final double demographicAtCoalPoint = demographicFunction.getPopSize(finishTime);

                    if (duration == 0.0 || demographicAtCoalPoint * (intervalArea / duration) >= threshold) {
                        logL -= Math.log(demographicAtCoalPoint);
                    } else {
                        return Double.NEGATIVE_INFINITY;
                    }

                } else {
                    double numerator = Math.exp(-kChoose2 * intervalArea) - Math.exp(-kChoose2 * normalisationArea);
                    logL += Math.log(numerator);

                }

                // normalisation

                double normExp = Math.exp(-kChoose2 * normalisationArea);

                double logDenominator;

                // the denominator has an irritating tendency to round to zero

                if (normExp != 1) {
                    logDenominator = Math.log1p(-normExp);
                } else {
                    logDenominator = handleDenominatorUnderflow(-kChoose2 * normalisationArea);
                }


                logL -= logDenominator;

            }

            startTime = finishTime;

        }



        return logL;
    }


    private static double handleDenominatorUnderflow(double input){
        BigDecimal bigDec = new BigDecimal(input);
        BigDecimal expBigDec = BigDecimalUtils.exp(bigDec, bigDec.scale());
        BigDecimal one = new BigDecimal(1.0);
        BigDecimal oneMinusExpBigDec = one.subtract(expBigDec);
        BigDecimal logOneMinusExpBigDec = BigDecimalUtils.ln(oneMinusExpBigDec, oneMinusExpBigDec.scale());
        return logOneMinusExpBigDec.doubleValue();
    }
}


