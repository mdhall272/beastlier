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
import beast.core.Input;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.StandardData;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public class Outbreak extends BEASTObject implements DataType {

    public Input<Boolean> hasGeographyInput = new Input<>("hasGeography", "True if hosts have geographical locations",
            false, Input.Validate.OPTIONAL);
    public Input<List<ClinicalCase>> casesInput = new Input<>("cases", "A list if all clinical cases (hosts) in the " +
            "outbreak; can include never-infected susceptibles");

    private Boolean hasGeography;
    private List<ClinicalCase> cases;
    private List<ClinicalCase> everInfectedCases;

    @Override
    public void initAndValidate() {
        cases = casesInput.get();
        hasGeography = hasGeographyInput.get();
        everInfectedCases = new ArrayList<>();

        for(ClinicalCase aCase : cases){
            if(aCase.wasEverInfected()){
                everInfectedCases.add(aCase);
            }
            if(hasGeography && !(aCase instanceof GeographicallyLocatedClinicalCase)){
                throw new IllegalArgumentException("Clinical case "+aCase.getID()+" is not geographically located" +
                        " but a geographical component is specified by outbreak "+getID());
            }
        }

    }

    @Override
    public int getStateCount() {
        return everInfectedCases.size();
    }

    @Override
    public List<Integer> string2state(String sequence) {
        return null;
    }

    @Override
    public String state2string(List<Integer> states) {
        return null;
    }

    @Override
    public String state2string(int[] states) {
        return null;
    }

    @Override
    public boolean[] getStateSet(int state) {
        return new boolean[0];
    }

    @Override
    public int[] getStatesForCode(int state) {
        return new int[0];
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
        return null;
    }

    @Override
    public char getChar(int state) {
        return 0;
    }

    @Override
    public String getCode(int state) {
        return null;
    }
}
