/**
 * @file ScoredPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.ScoredPeptideHit
 * @author Tao Xu
 * @date $Date: 2007/04/12 22:49:56 $
 */



package edu.scripps.pms.topdown;

import edu.scripps.pms.util.seq.Fasta;
import java.util.*;

public class ScoredPeptideHit implements Comparable<ScoredPeptideHit> {
    // Each ScoredPeptideHit may be associated with multiple PeptideHit
    private double xcorr;
    private double pScore; // probability score
    private int xcorrRank;
    private int pScoreRank;
    private LinkedList<PeptideHit> hits = new LinkedList<PeptideHit>();
    private String sequence;
    
    // 0 for sort by binomial probability,
    // 1 for sort by xcorr 
    private int additionalEstimate;


    public int [] getTheorMasses() {
        return hits.get(0).getTheorMasses();
    }
    public double getTheorMass() {
        return hits.get(0).getTheorMass();
    }
    public ScoredPeptideHit(String seq, int addE) {
        sequence = seq;
        this.additionalEstimate = addE;
    }
    public boolean isModified() {
        return hits.get(0).isModified();
    }
    public double getNTermMassShift() {
        return hits.get(0).getNTermMassShift();
    }
    public double getCTermMassShift() {
        return hits.get(0).getCTermMassShift();
    }
    public String getSequence() {
        return sequence;
    }
    public String getExtendedSequence() {
        return hits.get(0).getExtendedSequence();
    }
    public int getSecondaryRank() {
        switch(additionalEstimate) {
            case 1: return pScoreRank;
            case 0: return xcorrRank; 
            default: return xcorrRank;
        }    
    }
    public double getSecondaryScore() {
        switch(additionalEstimate) {
            case 1: return pScore;
            case 0: return xcorr; 
            default: return xcorr;
        }    
    }
    public int getPrimaryRank() {
        switch(additionalEstimate) {
            case 0: return pScoreRank;
            case 1: return xcorrRank; 
            default: return pScoreRank;
        }    
    }
    public double getPrimaryScore() {
        switch(additionalEstimate) {
            case 0: return pScore;
            case 1: return xcorr; 
            default: return pScore;
        }    
    }
    public void setPrimaryRank(int r) {
        
        switch(additionalEstimate) {
            case 0: pScoreRank = r;
            case 1: xcorrRank = r; 
            default: pScoreRank = r;
        }    
    }
    public void setXCorr(double x) {
        xcorr = x;
    }
    public void setPScore(double x) {
        pScore = x;
    }
    public void setPScoreRank(int r) {
        pScoreRank = r;
    }
    public void setXCorrRank(int r) {
        xcorrRank = r;
    }
    public double getPScore() {
        return pScore;
    }
    public double getXCorr() {
        return xcorr;
    }
    public int getNumPeaks() {
        return hits.get(0).getNumPeaks();
    }
    public int getNumPeaksMatched() {
        return hits.get(0).getNumPeaksMatched();
    } 
    public void addPeptideHit(PeptideHit p) {
        hits.add(p);
    } 
    public List<PeptideHit> getPeptideHits() {
        return hits;
    }

    // sort by xcorr 
    public int compareTo(ScoredPeptideHit s) {
        //double f = s.xcorr- this.xcorr;  // difference
        double f = s.getNumPeaksMatched() - this.getNumPeaksMatched();  // difference

        if (f > 0) {
            return 1;
        } else if (f < 0) {
            return -1;
        } else {
            return 0;
        }

    }
}
