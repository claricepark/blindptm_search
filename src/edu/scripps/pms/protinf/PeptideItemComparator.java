/**
 * @file PeptideItemComparater.java
 * This is the source file for edu.scripps.pms.protinf.PeptideItemComparater
 * @author Tao Xu
 * @date $Date
 */

package edu.scripps.pms.protinf;


import java.util.Comparator;
//import java.util.List;

public class PeptideItemComparator implements Comparator {

    // how to comapare two PeptideItemComparaters
    private int compareMode = 1; // 1 sort by number of PeptideItems, 2 sort by sumZScore, 3 by ProteinScore

    public PeptideItemComparator(int compareMode) {
        this.compareMode = compareMode;
    }   
    public void setCompareMode(int compareMode) {
        this.compareMode = compareMode;
    }


    public int compare(Object point1, Object point2) {

        PeptideItem p1 = (PeptideItem) point1;
        PeptideItem p2 = (PeptideItem) point2;
        double f = 0;  // difference
         

        switch(compareMode) {
            case 8 : f = p1.getAvgZScore() - p2.getAvgZScore(); break; // sort by average Z Score
            case 7 : f = p1.getConfidenceProduct() - p2.getConfidenceProduct(); break; // sort by confidence product
            case 6 : f = p1.getConfidenceSum() - p2.getConfidenceSum(); break; // sort by confidence sum 
            case 5 : f = p1.getBestXCorr() - p2.getBestXCorr(); break;
            //case 4 : f = p1.getOccurance() - p2.getOccurance(); break; // sort by number of peptides 
            case 4 : f = p1.getOccurrence() - p2.getOccurrence(); break; // sort by number of peptides 
            case 3 : f = p1.getSumXCorr() - p2.getSumXCorr(); break; // sort by sumXScore
            case 2 : f = p1.getSumZScore() - p2.getSumZScore(); break; // sort by sumZScore
            case 1 : f = p1.getBestZScore() - p2.getBestZScore(); break;
            default : f = p1.getBestZScore() - p2.getBestZScore(); // sort by protein length
        }


        
        if (f > 0) {
            return -1;
        } else if (f < 0) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * Return true if both the intensities and xValues are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {
        PeptideItemComparator c = (PeptideItemComparator)o;
        if(o == null) {
            return false;
        }
        return this.compareMode == c.compareMode;
    }
}
