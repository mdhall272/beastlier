/*
* File BetweenHostModel.java
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
import beast.evolution.tree.EpidemiologicalPartitionedTree;
import beast.evolution.tree.TreeDistribution;
import beastlier.outbreak.ClinicalCase;
import beastlier.outbreak.Outbreak;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public abstract class BetweenHostModel extends TreeDistribution {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The collection of clinical cases");

    protected EpidemiologicalPartitionedTree tree;
    protected Outbreak outbreak;
    protected boolean hasLatentPeriods;

    protected ArrayList<TreeEvent> sortedTreeEvents;
    protected ArrayList<TreeEvent> storedSortedTreeEvents;

    public enum EventType{
        INFECTION,
        INFECTIOUSNESS,
        END
    }

    public void initAndValidate(){

        if(!(treeInput.get() instanceof EpidemiologicalPartitionedTree)){
            throw new IllegalArgumentException("Trees given to the between-host model must have node partitions and" +
                    " an outbreak");
        }

        tree = (EpidemiologicalPartitionedTree) treeInput.get();
        outbreak = outbreakInput.get();
        sortEvents();
        storedSortedTreeEvents = new ArrayList<>(sortedTreeEvents);
    }

    protected double getInfectionTime(ClinicalCase aCase){
        return tree.getInfectionTime(aCase);
    }
    protected abstract double getInfectiousTime(ClinicalCase aCase);

    protected class TreeEvent{

        private EventType type;
        private double time;
        private ClinicalCase aCase;
        private ClinicalCase infectorCase;

        private TreeEvent(EventType type, double time, ClinicalCase aCase){
            this.type = type;
            this.time = time;
            this.aCase = aCase;
            this.infectorCase = null;
        }

        private TreeEvent(double time, ClinicalCase aCase, ClinicalCase infectorCase){
            this.type = EventType.INFECTION;
            this.time = time;
            this.aCase = aCase;
            this.infectorCase = infectorCase;
        }

        public double getTime(){
            return time;
        }

        public EventType getType(){
            return type;
        }

        public ClinicalCase getCase(){
            return aCase;
        }

        public ClinicalCase getInfector(){
            return infectorCase;
        }

    }


    protected void sortEvents(){
        ArrayList<TreeEvent> out = new ArrayList<>();
        for(ClinicalCase aCase : outbreak.getCases()){

            double infectionTime = Double.POSITIVE_INFINITY;

            if(aCase.wasEverInfected()){
                infectionTime = tree.getInfectionTime(aCase);
            }

            out.add(new TreeEvent(infectionTime, aCase, tree.getInfector(aCase)));

            if(aCase.wasEverInfected()) {

                double endTime = aCase.getEndTime();

                out.add(new TreeEvent(EventType.END, endTime, aCase));

                if (hasLatentPeriods) {
                    double infectiousnessTime = getInfectiousTime(aCase);
                    out.add(new TreeEvent(EventType.INFECTIOUSNESS, infectiousnessTime, aCase));

                }
            }
        }

        Collections.sort(out, new EventComparator());

        sortedTreeEvents = out;

    }

    private class EventComparator implements Comparator<TreeEvent> {
        public int compare(TreeEvent treeEvent1, TreeEvent treeEvent2) {
            return Double.compare(treeEvent1.getTime(),
                    treeEvent2.getTime());
        }
    }


}
