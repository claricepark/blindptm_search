/**
 * @file ProteinHit.java
 * This is the source file for edu.scripps.pms.util.seq.ProteinHit
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.blindssptm;



import edu.scripps.pms.util.seq.Fasta;

public class ProteinHit implements Comparable<ProteinHit> {

    // description line of this ProteinHit sequence
    // the sequence string of this ProteinHit
    
    //private int chargeState;
    private int numPeaksMatched = 0;
    private int numTheroticPeaks = 0;
    private double probability;
    private Fasta protein;

    public ProteinHit(Fasta protein) {
        this.protein = protein;
    }  
    
    public boolean equals(ProteinHit p) {
        
        if(p != null) {
            return getSequence().equals(p.getSequence()); 
        }
        return false;
    }
    public int compareTo(ProteinHit p) {

//        return p.numPeaksMatched - this.numPeaksMatched; 
        if(p.probability < this.probability) {
            return 1;
        } else if(p.probability > this.probability) {
            return -1;
        } else { return 0; }
       
    }
    public void setNumPeaksMatched(int n, int t, double ptrue) {
        numPeaksMatched = n;
        numTheroticPeaks = t;
        probability = DistributionCalculator.calcBinomialSum(ptrue, t, n);
    }
    public Fasta getProtein() {
        return protein; 
    }
    public double getProbability() {
        return probability;
    }
    public String getSequence() {
        return protein.getSequence(); 
    }
    public String getDefline() {
        return protein.getDefline();
    }

    public int getNumPeaksMatched() {
        return numPeaksMatched;
    }
    public int getLength() {
        return protein.getLength();
    }
    public String getAccession() {
        return protein.getAccession();
    }
}
