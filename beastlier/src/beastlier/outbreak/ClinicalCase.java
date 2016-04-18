/*
* File ClinicalCase.java
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
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

@Description("Base class for epidemiological information about a case in an outbreak")
public class ClinicalCase extends BEASTObject {

    public Input<TaxonSet> taxaInput = new Input<>("taxa", "The taxa taken from this host");
    public Input<RealParameter> endOfInfectiousTimeInput = new Input<>("endOfInfectiousTime", "The time at which" +
            " this host ceased to be infectious", null, Input.Validate.OPTIONAL);
    public Input<Boolean> wasEverInfectedInput = new Input<>("wasEverInfected", "Was this host ever infected?", true,
            Input.Validate.OPTIONAL);
    public Input<Double> indexPriorWeightInput = new Input<>("indexPriorWeight", "The prior weight placed on this " +
            "host being the index infection of the outbreak", 1.0, Input.Validate.OPTIONAL);

    private TaxonSet taxa;
    private RealParameter endOfInfectiousTime;
    private Boolean wasEverInfected;
    private double indexPriorWeight;

    public boolean wasEverInfected(){
        return wasEverInfected;
    }

    public double getEndTime(){
        return endOfInfectiousTime.getValue();
    }

    public double getIndexPriorWeight(){
        return indexPriorWeight;
    }

    @Override
    public void initAndValidate() {
        taxa = taxaInput.get();
        endOfInfectiousTime = endOfInfectiousTimeInput.get();
        wasEverInfected = wasEverInfectedInput.get();
        indexPriorWeight = indexPriorWeightInput.get();

        if(wasEverInfected && endOfInfectiousTime==null){
            throw new IllegalArgumentException("Case "+this.getID()+" not given as never-infected but has no " +
                    "end of infectiousness parameter" );
        }

        if(!wasEverInfected && endOfInfectiousTime!=null){
            throw new IllegalArgumentException("Case "+this.getID()+" given as never-infected but has specified" +
                    " end of infectiousness parameter" );
        }
    }
}
