/**
 * @file ProteinGroupComparater.java
 * This is the source file for edu.scripps.pms.protinf.ProteinGroupComparater
 * @author Tao Xu
 * @date $Date
 */

package edu.scripps.pms.protinf;


import java.util.Comparator;
//import java.util.List;

public class ProteinGroupComparator implements Comparator {

    // how to comapare two ProteinGroupComparaters
    private int compareMode = 1; // 1 sort by number of PeptideItems, 2 sort by sumZScore, 3 by ProteinScore

    public ProteinGroupComparator(int compareMode) {
        this.compareMode = compareMode;
    }   
    public void setCompareMode(int compareMode) {
        this.compareMode = compareMode;
    }
    private double sortByNumTrypticPeptides (ProteinGroup p1, ProteinGroup p2) {
        double diff = p1.getIdentifiedTrypticPeptideFraction() - p2.getIdentifiedTrypticPeptideFraction();
        if(diff == 0) {
            diff = p1.getRepresentative().getAverageZScore() - p2.getRepresentative().getAverageZScore();
        }
        return diff;
    }

    public int compare(Object point1, Object point2) {

        ProteinGroup p1 = (ProteinGroup) point1;
        ProteinGroup p2 = (ProteinGroup) point2;
        double f = 0;  // difference
         

        switch(compareMode) {
            case 8 : f = p1.getIdentifiedTrypticPeptideFraction() - p2.getIdentifiedTrypticPeptideFraction(); break; // sort by numIdentifiedTyprticPeptides/numTrypticPeptide then average Z Score, then protein length
            case 7 : f = p1.getConfidenceProduct() - p2.getConfidenceProduct(); break; // sort by confidence product
            case 6 : f = p1.getConfidenceSum() - p2.getConfidenceSum(); break; // sort by confidence sum 
            case 5 : f = p1.getRepresentative().getAverageZScore() - p2.getRepresentative().getAverageZScore(); break; // sort by average Z Score 
            case 4 : f = p1.getRepresentative().getNumPeptides() - p2.getRepresentative().getNumPeptides(); break; // sort by number of peptides 
            case 3 : f = p1.getRepresentative().getProteinScore() - p2.getRepresentative().getProteinScore(); break; // sort by sumZScore
            case 2 : f = p1.getRepresentative().getSumZScore() - p2.getRepresentative().getSumZScore(); break; // sort by sumZScore
            case 1 : f = p1.getRepresentative().getFasta().getLength() - p2.getRepresentative().getFasta().getLength(); break;

            default : f = p1.getRepresentative().getFasta().getLength() - p2.getRepresentative().getFasta().getLength(); break; // sort by protein length
        }


        
        if (f > 0) {
            return -1;
        } else if (f < 0) {
            return 1;
        } else {
            if(p1.getNumPeptideItems() > p2.getNumPeptideItems()) {
                return -1;
            }  else if(p1.getNumPeptideItems() < p2.getNumPeptideItems()) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    /**
     * Return true if both the intensities and xValues are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {
        ProteinGroupComparator c = (ProteinGroupComparator)o;
        if(o == null) {
            return false;
        }
        return this.compareMode == c.compareMode;
    }
}
