/*
* File ToySI.java
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

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;
import beastlier.durations.DurationCategory;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.GeographicallyLocatedClinicalCase;
import org.apache.commons.math.FunctionEvaluationException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class ToySI extends BetweenHostModel {

    public Input<RealParameter> transmissionRateInput = new Input<>("transmissionRate", "The transmission rate");

    public Input<ParametricDistribution> initialInfectionTimePriorInput = new Input<>("initialInfectionTimePrior",
            "The prior distribution for the time of the index infection", null, Input.Validate.OPTIONAL);


    RealParameter transmissionRate;
    ParametricDistribution initialInfectionTimePrior;
    private HashMap<ClinicalCase, Double> indexCasePrior;

    public void initAndValidate(){

        hasLatentPeriods = false;


        transmissionRate = transmissionRateInput.get();
        initialInfectionTimePrior = initialInfectionTimePriorInput.get();

        double totalWeights = 0;
        super.initAndValidate();

        indexCasePrior = new HashMap<>();

        for(ClinicalCase aCase : outbreak.getCases()){
            totalWeights += aCase.getIndexPriorWeight();
        }

        for(ClinicalCase aCase : outbreak.getCases()){
            indexCasePrior.put(aCase, aCase.getIndexPriorWeight()/totalWeights);
        }
    }

    @Override
    public double evaluateLogP() {
        double transLogProb = 0;

        double rate = transmissionRate.getValue();

        ArrayList<ClinicalCase> previouslyInfectious = new ArrayList<>();

        double currentEventTime;
        boolean first = true;

        for (TreeEvent event : sortedTreeEvents) {
            currentEventTime = event.getTime();

            ClinicalCase thisCase = event.getCase();

            if (event.getType() == EventType.INFECTION) {
                if (first) {
                    // index infection

                    if (indexCasePrior != null) {
                        transLogProb += Math.log(indexCasePrior.get(thisCase));
                    }
                    if (initialInfectionTimePrior != null) {
                        transLogProb += initialInfectionTimePrior.logDensity(currentEventTime);
                    }
                    if (!hasLatentPeriods) {
                        previouslyInfectious.add(thisCase);
                    }

                    first = false;

                } else {

                    ClinicalCase infector = event.getInfector();

                    if (thisCase.wasEverInfected()) {
                        if (previouslyInfectious.contains(thisCase)) {
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (event.getTime() > thisCase.getEndTime()) {
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (infector.getEndTime() < event.getTime()) {
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (getInfectiousTime(infector) > event.getTime()) {
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (!previouslyInfectious.contains(infector)) {
                            throw new RuntimeException("Infector not previously infected");
                        }
                    }

                    // no other previously infectious case has infected this case...

                    for (ClinicalCase nonInfector : previouslyInfectious) {

                        double timeDuringWhichNoInfection;
                        if (nonInfector.getEndTime() < event.getTime()) {
                            timeDuringWhichNoInfection = nonInfector.getEndTime() - getInfectiousTime(nonInfector);
                        } else {
                            timeDuringWhichNoInfection = event.getTime() - getInfectiousTime(nonInfector);
                        }

                        if (timeDuringWhichNoInfection < 0) {
                            throw new RuntimeException("negative time");
                        }

                        double transRate = rate;

                        transLogProb += -transRate * timeDuringWhichNoInfection;
                    }

                    // ...until the end

                    if (thisCase.wasEverInfected()) {
                        double transRate = rate;
                        transLogProb += Math.log(transRate);
                    }
                    if (!hasLatentPeriods) {
                        previouslyInfectious.add(thisCase);
                    }
                }
            }
        }


        // just reject states where these round to +INF

        if(transLogProb == Double.POSITIVE_INFINITY){
            System.out.println("TransLogProb +INF");
            return Double.NEGATIVE_INFINITY;
        }


        logP = transLogProb;

        return logP;
    }


    @Override
    protected double getInfectiousTime(ClinicalCase aCase) {
        return getInfectionTime(aCase);
    }
}
