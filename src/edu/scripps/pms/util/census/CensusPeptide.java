package edu.scripps.pms.util.census;

/**
 * @author Tao Xu    
 * @version $Id
 */
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

public class CensusPeptide {

    private String peptideline;
    private String sequence;
    private double ratio;
    private double regressionFactor; 
    private boolean isModified = false;
 
    public CensusPeptide(String line) {
        peptideline = line;
        String [] arr = line.split("\t");
        sequence = arr[2];
        ratio = Double.parseDouble(arr[3]);
        regressionFactor = Double.parseDouble(arr[4]);

        if(sequence.indexOf("*") > -1 || sequence.indexOf("#") > -1 || 
               sequence.indexOf("@") > -1 || sequence.indexOf("(") > -1) {
            isModified = true;
        }
    }
   
    public double getRatio() {
        return ratio;
    }

    public boolean isModifiedPeptide() {
        return isModified;
    }    
}
