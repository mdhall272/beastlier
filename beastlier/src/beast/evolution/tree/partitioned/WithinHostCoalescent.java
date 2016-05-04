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

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.PartitionedTreeNode;
import beast.evolution.tree.coalescent.*;
import beast.math.Binomial;
import beastlier.outbreak.ClinicalCase;
import beast.util.BigDecimalUtils;

import java.io.PrintStream;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("A parameteric within-host coalescent function shared by every clinical case")
public class WithinHostCoalescent extends WithinHostModel {

    public Input<PopulationFunction> functionInput = new Input<>("populationFunction", "The within-host coalescent " +
            "process");

    private PopulationFunction popFunction;
    protected static double tolerance = 1E-10;
    public boolean[] recalculateTreeletLogP;
    private double[] individualWHProbabilities;
    private double[] storedIndividualWHProbabilities;


    public void initAndValidate(){
        super.initAndValidate();
        popFunction = functionInput.get();
        individualWHProbabilities = new double[tree.getElementList().size()];
        recalculateTreeletLogP = new boolean[tree.getElementList().size()];
        Arrays.fill(recalculateTreeletLogP, true);
    }


    public double calculateLogP() {

        RealParameter q = tree.getQ();
        for (int i = 0; i < q.getDimension(); i++) {
            if (q.isDirty(i)) {
                treeletsRequiringExtraction[i] = true;
                recalculateTreeletLogP[i] = true;
                int parentPartitionNumber = tree.getAncestorPartitionElement(i);
                if (parentPartitionNumber != -1) {
                    treeletsRequiringExtraction[parentPartitionNumber] = true;
                    recalculateTreeletLogP[parentPartitionNumber] = true;
                }
            }
        }

        for (Node node : tree.getNodesAsArray()) {
            PartitionedTreeNode castNode = (PartitionedTreeNode) node;
            if (castNode.isPartitionDirty()) {
                treeletsRequiringExtraction[castNode.getPartitionElementNumber()] = true;
                recalculateTreeletLogP[castNode.getPartitionElementNumber()] = true;
            }
        }

        // if the population function has changed, then all treelets need probabilities recalculated but (unless
        // something else has changed) no treelets actually need re-extracting

        if(((CalculationNode)popFunction).isDirtyCalculation()){
            Arrays.fill(recalculateTreeletLogP, true);
        }

        explodeTree();

        logP = 0;

        for (int i=0; i<tree.getNElements(); i++) {

            if(recalculateTreeletLogP[i]) {

                String caseName = tree.getElementString(i);
                ClinicalCase aCase = outbreak.getCaseByID(caseName);

                Treelet treelet = elementsAsTrees.get(aCase);

                if (treelet.getLeafNodeCount() > 1) {
                    TreeIntervals intervals = new TreeIntervals(treelet);

                    double individualLogP = calculateTreeletLogLikelihood(intervals, popFunction, 0,
                            treelet.getZeroHeight());

                    individualWHProbabilities[i] = individualLogP;

                    logP += individualLogP;

                } else {
                    individualWHProbabilities[i] = 0.0;

                    logP += 0.0;
                }
            } else {
                logP += individualWHProbabilities[i];
            }



        }

        return logP;
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

    @Override
    public void store() {
        storedElementsAsTrees = new HashMap<>(elementsAsTrees);
        storedIndividualWHProbabilities = individualWHProbabilities.clone();
        Arrays.fill(treeletsRequiringExtraction, false);
        Arrays.fill(recalculateTreeletLogP, false);
        super.store();
    }

    @Override
    public void restore() {
        elementsAsTrees = storedElementsAsTrees;
        individualWHProbabilities = storedIndividualWHProbabilities;
        Arrays.fill(treeletsRequiringExtraction, false);
        Arrays.fill(recalculateTreeletLogP, false);
        super.restore();
    }

    public void init(final PrintStream out){
        for(int i=0; i<tree.getElementList().size(); i++){
            out.print("whc_logP_"+tree.getElementList().get(i) + "\t");
        }
        out.print("whc_logP_total" + "\t");
    }

    public void log(final int sample, final PrintStream out){
        for(int i=0; i<tree.getElementList().size(); i++){
            out.print(individualWHProbabilities[i] + "\t");
        }
        out.print(logP + "\t");
    }

    private static double handleDenominatorUnderflow(double input){
        BigDecimal bigDec = new BigDecimal(input);
        BigDecimal expBigDec = BigDecimalUtils.exp(bigDec, bigDec.scale());
        BigDecimal one = new BigDecimal(1.0);
        BigDecimal oneMinusExpBigDec = one.subtract(expBigDec);
        BigDecimal logOneMinusExpBigDec = BigDecimalUtils.ln(oneMinusExpBigDec, oneMinusExpBigDec.scale());
        return logOneMinusExpBigDec.doubleValue();
    }

    @Override
    protected boolean requiresRecalculation() {

        return ((CalculationNode) functionInput.get()).isDirtyCalculation()
                || super.requiresRecalculation();
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    public List<String> getArguments(){
        return Collections.singletonList(treeInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    public List<String> getConditions(){
        return functionInput.get().getParameterIds();
    }

}


