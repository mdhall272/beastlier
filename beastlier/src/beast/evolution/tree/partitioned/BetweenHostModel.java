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
import beastlier.util.PartitionedTreeLogger;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Matthew Hall <mdhall@ic.ac.uk>
 */

public abstract class BetweenHostModel extends TreeDistribution {

    public Input<Outbreak> outbreakInput = new Input<>("outbreak", "The collection of clinical cases");

    protected EpidemiologicalPartitionedTree tree;
    protected Outbreak outbreak;
    protected boolean hasLatentPeriods;

    //Dirty here means model parameters have changed but the tree and its timings have not so the order of events does
    // not need recalculating. Filthy means they has so it does;

    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;
    protected int typeOfDirt = 2;

    protected ArrayList<TreeEvent> sortedTreeEvents;
    protected ArrayList<TreeEvent> storedSortedTreeEvents;
    protected AbstractMap<ClinicalCase, Double> infectionTimesMap;
    protected AbstractMap<ClinicalCase, Double> storedInfectionTimesMap;

    public enum EventType{
        INFECTION,
        INFECTIOUSNESS,
        END
    }

    // Subclasses _must_ do this first.

    public double calculateLogP(){

        if(typeOfDirt == IS_FILTHY) {
            infectionTimesMap = null;
            sortEvents();
        }

        return evaluateLogP();
    }

    //Abstract method for doing the actual calculations

    public abstract double evaluateLogP();

    public void initAndValidate(){

        if(!(treeInput.get() instanceof EpidemiologicalPartitionedTree)){
            throw new IllegalArgumentException("Trees given to the between-host model must have node partitions and" +
                    " an outbreak");
        }

        tree = (EpidemiologicalPartitionedTree) treeInput.get();
        outbreak = outbreakInput.get();
        sortEvents();
        storedSortedTreeEvents = new ArrayList<>(sortedTreeEvents);

        infectionTimesMap = null;
        storedInfectionTimesMap = null;
    }

    protected double getInfectionTime(ClinicalCase aCase){
        if(infectionTimesMap==null){
            getInfectionTimesFromTree();
        }
        return infectionTimesMap.get(aCase);
    }
    protected abstract double getInfectiousTime(ClinicalCase aCase);

    private void getInfectionTimesFromTree(){
        infectionTimesMap = new HashMap<>();
        for(ClinicalCase aCase : outbreak.getCases()){
            infectionTimesMap.put(aCase, tree.getInfectionTime(aCase));
        }
    }

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
//
//        try {
//            PrintStream firstSteam = new PrintStream("ttg.nex");
//
//            PartitionedTreeLogger.debugLog(tree, 0, false, firstSteam);
//        } catch (FileNotFoundException e){
//            e.printStackTrace();
//        }

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
            return Double.compare(treeEvent1.getTime(), treeEvent2.getTime());
        }
    }

    public void store() {
        storedInfectionTimesMap = infectionTimesMap == null ? null : new HashMap<>(infectionTimesMap);
        storedSortedTreeEvents = new ArrayList<>(sortedTreeEvents);
        typeOfDirt = IS_CLEAN;
        super.store();
    }


    public void restore() {
        infectionTimesMap = storedInfectionTimesMap;
        sortedTreeEvents = storedSortedTreeEvents;
        typeOfDirt = IS_CLEAN;
        super.restore();
    }

    // stop BEAUTi from giving within host models as combo box options

    public boolean canSetXx(Object o){
        return !(o instanceof WithinHostModel);
    }

    @Override
    protected boolean requiresRecalculation() {
        return treeInput.get().somethingIsDirty();
    }

}
