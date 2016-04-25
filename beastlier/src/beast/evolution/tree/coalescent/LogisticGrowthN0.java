/*
 * LogisticGrowthN0.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

/*
 * LogisticGrowthN0.java
 *
 * Daniel Wilson 4th October 2011
 *
 */

package beast.evolution.tree.coalescent;


import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.LambertW;

import java.util.ArrayList;
import java.util.List;

/**
 * This class models logistic growth. Ported directly from BEAST 1.
 *
 * @author Daniel Wilson
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Matthew Hall
 * @version $Id: LogisticGrowth.java,v 1.15 2008/03/21 20:25:56 rambaut Exp $
 */
public class LogisticGrowthN0 extends PopulationFunction.Abstract {

    final public Input<RealParameter> t50ParameterInput = new Input<>("t50",
            "Time of half");
    final public Input<RealParameter> growthRateParameterInput = new Input<>("growthRate",
            "growth rate is the exponent of the exponential growth");
    final public Input<RealParameter> popSizeParameterInput = new Input<>("popSize",
            "present-day population size (defaults to 1.0). ");

    @Override
    public void initAndValidate() {
        if (popSizeParameterInput.get() != null) {
            popSizeParameterInput.get().setBounds(
                    Math.max(0.0, popSizeParameterInput.get().getLower()),
                    popSizeParameterInput.get().getUpper());
        }
//        if (growthRateParameter.get() != null) {
//            growthRateParameter.get().setBounds(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//        }
    }

    public double getN0() {
        return popSizeParameterInput.get().getValue();
    }
    public final double getGrowthRate() {
        return growthRateParameterInput.get().getValue();
    }
    public final double getT50() {return t50ParameterInput.get().getValue(); }

    // Implementation of abstract methods

    @Override
    public List<String> getParameterIds() {
        ArrayList<String> out = new ArrayList<>();
        out.add(popSizeParameterInput.get().getID());
        out.add(growthRateParameterInput.get().getID());
        out.add(t50ParameterInput.get().getID());
        return out;
    }

    /**
     * Gets the value of the demographic function N(t) at time t.
     *
     * @param t the time
     * @return the value of the demographic function N(t) at time t.
     */
    public double getPopSize(double t) {
        double N0 = getN0();
        double r = getGrowthRate();
        double T50 = getT50();

		return N0 * (1 + Math.exp(-r * T50)) / (1 + Math.exp(-r * (T50-t)));
    }

    public double getLogDemographic(double t) {
		return Math.log(getPopSize(t));
	}

    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {
        double N0 = getN0();
        double r = getGrowthRate();
        double T50 = getT50();
		double exp_rT50 = Math.exp(-r*T50);

		return (t + exp_rT50 * (Math.exp(r * t) - 1)/r) / (N0 * (1 + exp_rT50));
    }

    /**
     * Returns the inverse function of getIntensity
     *
     * If exp(-qt) = a(t-k) then t = k + (1/q) * W(q*exp(-q*k)/a) where W(x) is the Lambert W function.
     *
     * for our purposes:
     *
     * q = -r
     * a = (q/exp(q*T50))
     * k = N0*(1+exp(q*T50))*x - (1/q)exp(q*T50)
     *
     * For large x, W0(x) is approximately equal to ln(x) - ln(ln(x)); if q*exp(-q*k)/a rounds to infinity, we log it
     * and use this instead
     */
    public double getInverseIntensity(double x) {

        double q = -getGrowthRate();
        double T50 = getT50();
        double N0 = getN0();
        double a = (q/Math.exp(q*T50));
        double k = N0*(1+Math.exp(q*T50))*x - (1/q)*Math.exp(q*T50);

        double lambertInput = q*Math.exp(-q*k)/a;

        double lambertResult;

        if(lambertInput==Double.POSITIVE_INFINITY){

            //use the asymptote; note q/a = exp(q*T50)

            double logInput = q*T50-q*k;
            lambertResult = logInput - Math.log(logInput);

        } else {
            lambertResult = LambertW.branch0(lambertInput);
        }

        return k + (1/q)*lambertResult;

    }

    public double getIntegral(double start, double finish) {
		return getIntensity(finish) - getIntensity(start);
    }

    public int getNumArguments() {
        return 3;
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0:
                return "N0";
            case 1:
                return "r";
            case 2:
                return "t50";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN0();
            case 1:
                return getGrowthRate();
            case 2:
                return getT50();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }


}
