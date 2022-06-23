/**
 * @file ProteinItemComparater.java
 * This is the source file for edu.scripps.pms.protinf.ProteinItemComparater
 * @author Tao Xu
 * @date $Date
 */

package edu.scripps.pms.protinf;


import java.util.Comparator;
//import java.util.List;

public class ProteinItemComparator implements Comparator {

    // how to comapare two ProteinItemComparaters
    private int compareMode = 1; // 1 sort by protein length, 2 sort by sumZScore, 3 by ProteinScore

    public ProteinItemComparator(int compareMode) {
        this.compareMode = compareMode;
    }   
    public void setCompareMode(int compareMode) {
        this.compareMode = compareMode;
    }


    public int compare(Object point1, Object point2) {

        ProteinItem p1 = (ProteinItem) point1;
        ProteinItem p2 = (ProteinItem) point2;
        double f = 0;  // difference
         

        switch(compareMode) {
            case 4 : f = p1.getNumPeptides() - p2.getNumPeptides(); break; // sort by number of peptides 
            case 3 : f = p1.getProteinScore() - p2.getProteinScore(); break; // sort by sumZScore
            case 2 : f = p1.getSumZScore() - p2.getSumZScore(); break; // sort by sumZScore
            case 1 : f = p1.getFasta().getLength() - p2.getFasta().getLength(); break;
            default : f = p1.getFasta().getLength() - p2.getFasta().getLength(); // sort by protein length
        }


        
        if (f > 0) {
            return 1;
        } else if (f < 0) {
            return -1;
        } else {
            return 0;
        }
    }

    /**
     * Return true if both the intensities and xValues are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {
        ProteinItemComparator c = (ProteinItemComparator)o;
        if(o == null) {
            return false;
        }
        return this.compareMode == c.compareMode;
    }
}
