/*
* File Outbreak.java
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

package beastlier.outbreak;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("An outbreak is a collection of ClinicalCases from the same event")
public class Outbreak extends TaxonSet implements DataType {

    public Input<Boolean> hasGeographyInput = new Input<>("hasGeography", "True if hosts have geographical locations",
            false, Input.Validate.OPTIONAL);
    public Input<List<ClinicalCase>> casesInput = new Input<>("clinicalCase", "A clinical case (host) in the " +
            "outbreak; may be a susceptible that was never infected", new ArrayList<>());

    private Boolean hasGeography;
    private List<ClinicalCase> cases;

    @Override
    public void initAndValidate() {

        cases = casesInput.get();
        hasGeography = hasGeographyInput.get();

        //The taxonList is the set of cases that were ever infected

        taxonList = new ArrayList<>();

        for(ClinicalCase aCase : cases){
            if(aCase.wasEverInfected()){
                taxonList.add(aCase);
            }
            if(hasGeography && !(aCase instanceof GeographicallyLocatedClinicalCase)){
                throw new IllegalArgumentException("Clinical case "+aCase.getID()+" is not geographically located" +
                        " but a geographical component is specified by outbreak "+getID());
            }
        }
        taxaNames = new ArrayList<>();
        for (final Taxon taxon : taxonList) {
            taxaNames.add(taxon.getID());
        }

    }

    @Override
    public int getStateCount() {
        return taxonList.size();
    }

    // I assume that nobody in their right mind would want to overlay two or more partition models. A "sequence"
    // has just one position

    @Override
    public List<Integer> string2state(String sequence) {
        ArrayList<Integer> out = new ArrayList<>();
        for(ClinicalCase aCase : cases){
            if(aCase.getID().equals(sequence)){
                out.add(cases.indexOf(aCase));
            }
        }

        if(out.size() == 0){
            throw new RuntimeException("String "+sequence+" not found amongst case IDs");
        }

        return out;
    }

    @Override
    public String state2string(List<Integer> states) {
        if(states.size() > 1){
            throw new RuntimeException("Outbreak is a one position per tip DataType");
        }

        return cases.get(states.get(0)).getID();
    }

    @Override
    public String state2string(int[] states) {
        if(states.length > 1){
            throw new RuntimeException("Outbreak is a one position per tip DataType");
        }

        return cases.get(states[0]).getID();
    }

    @Override
    public boolean[] getStateSet(int state) {
        boolean[] out = new boolean[taxonList.size()];
        Arrays.fill(out, true);
        return out;
    }

    @Override
    public int[] getStatesForCode(int state) {
        int[] out = new int[1];
        out[0] = state;
        return out;
    }

    @Override
    public boolean isAmbiguousState(int state) {
        // no ambiguities allowed here

        return false;
    }

    @Override
    public boolean isStandard() {
        return false;
    }

    @Override
    public String getTypeDescription() {
        return "This DataType consists of all hosts/clinical cases in a known outbreak";
    }

    //These don't seem important

    @Override
    public char getChar(int state) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String getCode(int state) {
        throw new RuntimeException("Not implemented");
    }

    public ClinicalCase getCaseByID(String id){
        for(ClinicalCase aCase : cases){
            if(aCase.getID().equals(id)){
                return aCase;
            }
        }
        throw new RuntimeException("Looking for a clinical case that is not present in the outbreak");
    }

    public ClinicalCase getCaseByTaxon(String taxon){
        for(ClinicalCase aCase : cases){
            if(aCase.getTaxonSet().containsAll(Lists.newArrayList(taxon))){
                return aCase;
            }
        }
        throw new RuntimeException("Looking for a clinical case for a nonexistant taxon");
    }

    public List<ClinicalCase> getCases(){
        return cases;
    }

    public List<Taxon> getEverInfectedCases(){
        return taxonList;
    }

    public ClinicalCase getCase(int number){
        return cases.get(number);
    }

    public ClinicalCase getEverInfectedCase(int number){
        return (ClinicalCase)taxonList.get(number);
    }

    public int getCaseIndex(ClinicalCase aCase){
        return cases.indexOf(aCase);
    }

    public int getInfectedCaseIndex(ClinicalCase aCase){
        return taxonList.indexOf(aCase);
    }

    public boolean hasGeography(){
        return hasGeography;
    }

}
