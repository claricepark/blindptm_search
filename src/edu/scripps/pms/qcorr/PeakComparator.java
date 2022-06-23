/**
 * @file PeakComparater.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeakComparater
 * @author Tao Xu
 * @date $Date: 2007/02/27 01:31:02 $
 */



package edu.scripps.pms.qcorr;

import java.util.Comparator;
//import java.util.List;

public class PeakComparator implements Comparator {

    // how to comapare two PeakComparaters
    private static boolean compareByIntensity = false; 

    public PeakComparator(boolean compareByIntensity) {
        this.compareByIntensity = compareByIntensity;
    }   
    public void setCompareMode(boolean compareMode) {
        compareByIntensity = compareMode;
    }


    public int compare(Object peak1, Object peak2) {

        Peak p1 = (Peak) peak1;
        Peak p2 = (Peak) peak2;
        float f = 0;  // difference
         
        if (compareByIntensity) {
           f = p1.getIntensity() - p2.getIntensity();
        } else {
            f = p1.getM2z() - p2.getM2z();
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
     * Return true if both the intensities and m2zs are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {
        PeakComparator c = (PeakComparator)o;
        if(o == null) {
            return false;
        }
        return this.compareByIntensity == c.compareByIntensity;
    }
}
