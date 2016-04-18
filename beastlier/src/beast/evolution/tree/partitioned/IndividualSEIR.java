/*
* File IndividualSEIR.java
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

import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.CompoundValuable;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;
import beastlier.durations.DurationCategory;
import beastlier.durations.FixedValueDurationCategory;
import beastlier.geography.SpatialKernel;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.GeographicallyLocatedClinicalCase;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class IndividualSEIR extends BetweenHostModel {

    public Input<SpatialKernel> kernelInput = new Input<>("kernel", "The spatial kernel function; if null or absent, " +
            "geography will be ignored", null, Input.Validate.OPTIONAL);
    public Input<List<DurationCategory>> latentInput = new Input<>("latentInput", "One or more categories uniting " +
            "latent periods across clinical cases");
    public Input<List<DurationCategory>> infectiousInput = new Input<>("infectiousInput", "One or more categories " +
            "uniting infectious periods across clinical cases");
    public Input<RealParameter> baseTransmissionRateInput = new Input<>("baseTransmissionRate", "The unmodified " +
            "transmission rate");
    public Input<ParametricDistribution> initialInfectionTimePriorInput = new Input<>("initialInfectionTimePriod", "The prior" +
            " distribution for the time of the index infection");

    private SpatialKernel kernel;
    private RealParameter baseTransmissionRate;
    private ParametricDistribution initialInfectionTimePrior;
    private HashMap<ClinicalCase, Double> indexCasePrior;
    private boolean hasGeography;
    private List<DurationCategory> latentCategories;
    private List<DurationCategory> infectiousCategories;

    public void initAndValidate(){
        super.initAndValidate();
        kernel = kernelInput.get();
        baseTransmissionRate = baseTransmissionRateInput.get();
        initialInfectionTimePrior = initialInfectionTimePriorInput.get();

        latentCategories = latentInput.get();

        for(DurationCategory category : latentCategories){
            if(category.hasProbability()){
                throw new RuntimeException("Latent periods in this model are identical");
            }

            for(ClinicalCase aCase : category.getCases()){
                for(DurationCategory anotherCategory : latentCategories){
                    if(category != anotherCategory){
                        if(anotherCategory.getCases().contains(aCase)){
                            throw new RuntimeException("A clinical case is in two latent period categories");
                        }
                    }
                }
            }
        }

        infectiousCategories = infectiousInput.get();

        for(DurationCategory category : infectiousCategories){
            for(ClinicalCase aCase : category.getCases()){
                for(DurationCategory anotherCategory : infectiousCategories){
                    if(category != anotherCategory){
                        if(anotherCategory.getCases().contains(aCase)){
                            throw new RuntimeException("A clinical case is in two infectious period categories");
                        }
                    }
                }
            }
        }

        double totalWeights = 0;

        indexCasePrior = new HashMap<>();

        for(ClinicalCase aCase : outbreak.getCases()){
            totalWeights += aCase.getIndexPriorWeight();
        }

        for(ClinicalCase aCase : outbreak.getCases()){
            indexCasePrior.put(aCase, aCase.getIndexPriorWeight()/totalWeights);
        }

        if(kernel!=null && outbreak.hasGeography()){
            throw new IllegalArgumentException("Kernel specified but ClinicalCases have no attached geographical" +
                    " information");
        }

        hasGeography = kernel != null;

        hasLatentPeriods = true;

    }

    public double calculateLogP(){

        //todo likelihoodKnown works very differently to BEAST 1. For now, calculate everything every time

        double transLogProb = 0;

        if (sortedTreeEvents == null) {
            sortEvents();
        }

        double rate = baseTransmissionRate.getValue();

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
                        transLogProb += initialInfectionTimePrior.density(currentEventTime);
                    }

                    if (!hasLatentPeriods) {
                        previouslyInfectious.add(thisCase);
                    }

                    first = false;

                } else {

                    ClinicalCase infector = event.getInfector();

                    if(thisCase.wasEverInfected()) {

                        if (previouslyInfectious.contains(thisCase)){
//                                throw new RuntimeException(thisCase.getID() + " infected after it was infectious");
                            return Double.NEGATIVE_INFINITY;
                        }

                        if (event.getTime() > thisCase.getEndTime()){
//                                throw new RuntimeException(thisCase.getID() + " ceased to be infected before it was " +
//                                        "infected");
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (infector.getEndTime() < event.getTime()){
//                                throw new RuntimeException(thisCase.getID() + " infected by " + infector.getID() +
//                                        " after the latter ceased to be infectious");
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (getInfectiousTime(infector) > event.getTime()) {
//                                throw new RuntimeException(thisCase.getID() + " infected by "
//                                        + infector.getID() + " before the latter became infectious");
                            return Double.NEGATIVE_INFINITY;
                        }

                        if(!previouslyInfectious.contains(infector)){
                            throw new RuntimeException("Infector not previously infected");
                        }
                    }

                    // no other previously infectious case has infected this case...

                    for (ClinicalCase nonInfector : previouslyInfectious) {



                        double timeDuringWhichNoInfection;
                        if (nonInfector.getEndTime() < event.getTime()) {
                            timeDuringWhichNoInfection = nonInfector.getEndTime()
                                    - getInfectiousTime(nonInfector);
                        } else {
                            timeDuringWhichNoInfection = event.getTime()
                                    - getInfectiousTime(nonInfector);
                        }

                        if(timeDuringWhichNoInfection<0){
                            throw new RuntimeException("negative time");
                        }

                        double transRate = rate;
                        if (hasGeography) {
                            transRate *= kernel.getValue((GeographicallyLocatedClinicalCase)thisCase,
                                    (GeographicallyLocatedClinicalCase)nonInfector);
                        }

                        transLogProb += -transRate * timeDuringWhichNoInfection;


                    }

                    // ...until the end

                    if(thisCase.wasEverInfected()) {
                        double transRate = rate;


                        if (hasGeography) {
                            transRate *= kernel.getValue((GeographicallyLocatedClinicalCase)thisCase,
                                    (GeographicallyLocatedClinicalCase)infector);
                        }


                        transLogProb += Math.log(transRate);
                    }

                    if (!hasLatentPeriods) {
                        previouslyInfectious.add(thisCase);
                    }


                }


            } else if (event.getType() == EventType.INFECTIOUSNESS) {
                if (event.getTime() < Double.POSITIVE_INFINITY) {

                    if(event.getTime() > event.getCase().getEndTime()){
                        throw new RuntimeException(event.getCase().getID() + " noninfectious before" +
                                "infectious");
                    }

                    if (first) {
                        throw new RuntimeException("First event is not an infection");
                    }

                    previouslyInfectious.add(thisCase);
                }
            }
        }






        double periodsLogProb = 0;

        HashMap<String, ArrayList<Double>> infectiousPeriodsByCategory
                = new HashMap<String, ArrayList<Double>>();




        for(DurationCategory category : infectiousCategories){
            if(category.hasProbability()) {
                List<ClinicalCase> relevantCases = category.getCases();

                Double[] infectiousPeriods = new Double[relevantCases.size()];

                for (int i = 0; i < infectiousPeriods.length; i++) {
                    infectiousPeriods[i] = relevantCases.get(i).getEndTime() - getInfectiousTime(relevantCases.get(i));
                }

                RealParameter collectionOfValues = new RealParameter(infectiousPeriods);

                periodsLogProb += category.getLogProbability(collectionOfValues);
            }
        }

        // just reject states where these round to +INF

        if(transLogProb == Double.POSITIVE_INFINITY){
            System.out.println("TransLogProb +INF");
            return Double.NEGATIVE_INFINITY;
        }
        if(periodsLogProb == Double.POSITIVE_INFINITY){
            System.out.println("PeriodsLogProb +INF");
            return Double.NEGATIVE_INFINITY;
        }


        return periodsLogProb + transLogProb;

    }

    public double getInfectiousTime(ClinicalCase aCase){
        DurationCategory category = getLatentCategory(aCase);
        if(category.hasProbability()){
            throw new RuntimeException("Latent periods in this model are fixed");
        }

        FixedValueDurationCategory castCategory = (FixedValueDurationCategory)category;

        return getInfectionTime(aCase) + castCategory.getValue();

    }

    public DurationCategory getLatentCategory(ClinicalCase aCase){
        for(DurationCategory category : latentCategories){
            if(category.getCases().contains(aCase)){
                return category;
            }
        }
        throw new RuntimeException("Can't find a latent period category for case "+aCase.getID());
    }

    public DurationCategory getInfectiousCategory(ClinicalCase aCase){
        for(DurationCategory category : infectiousCategories){
            if(category.getCases().contains(aCase)){
                return category;
            }
        }
        throw new RuntimeException("Can't find an infectious period category for case "+aCase.getID());
    }


}
