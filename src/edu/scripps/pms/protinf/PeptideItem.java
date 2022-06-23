package edu.scripps.pms.protinf;

import java.util.*;
import edu.scripps.pms.util.dtaselect.Peptide;


/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class PeptideItem implements Comparable {

    private HashSet<ProteinItem> proteins = new HashSet();
    private ArrayList<ProteinGroup> protGroups = new ArrayList();
    private String seq;
    private int charge;
    private boolean isReverseHit = false;
    private Peptide [] peptides = null;
    
    private double sumXCorr = 0;
    private double sumZScore = 0;
    private double sumDeltaCn = 0; 
    private double bestZScore = 0;
    private double bestXCorr = 0;
    private int occurance = 0;
    private int totalSpectrumCount = 0;
    private int scoreType = 1; // 1 for bestZScore, 2 for sumZScore, 3 for sumXCorr, 4 for confidence
    
    private double xcorrConfidence = 0;
    private double zscoreConfidence = 0;
    private double confidence = 0; // confidence is calculated based on the xcorr confidence and confidence
   
    private double fdr = 100; // fdr of this peptide item in the dataset
   
    private ArrayList<Double> confidenceScores = new ArrayList();
    private int trypticStatus = 0;
 
    public void setTrypticStatus(int numend) {
        trypticStatus = numend;
    }
   
    public int getTrypticStatus() {
        return trypticStatus;
    }
    public String getPeptideInfo() {
        return seq + "\t" + charge + "\t" + trypticStatus + "\t" + getDecoyStatus() + "\t" + getPtmStatus() + "\t" + getSortingCode() + "\t" + getOccurrence() + "\t" + getBestPeptideScores() + "\t" + getTotalSpectrumCount() + "\t" + getConfidence() + "\t";
    } 
    public String getBestPeptideScores() {
        Peptide p = getBestPeptideHit();
        if(p == null) {
            return "-\t-\t-";
        } else {
            return p.getXCorr() + "\t" + p.getDeltCN() + "\t" + p.getSpScore() + "\t" + p.getDeltaMass() + "\t" + p.getFileNameScanCharge();
        }
    }
    public String getPtmStatus() {
        return seq.indexOf("(") != -1 || seq.indexOf("*") != -1 ||  seq.indexOf("@") != -1 || seq.indexOf("#") != -1? "TRUE" : "FALSE";
        
    }

    private String isDecoyHit() {
        boolean isdecoyhit = true;
System.out.println("Number of protein groups for " + seq + ": " + protGroups.size());
        for(Iterator<ProteinGroup> it = protGroups.iterator(); it.hasNext();) {
            ProteinGroup pg = it.next();
            isdecoyhit &= pg.isReverseHit();
            //if(seq.equals("KVAVEK")) {
System.out.println(seq + "\t" + pg.getRepresentativeAcc() + "\t" +  pg.isReverseHit());
            //} 
        }
        
        return isdecoyhit? "TRUE" : "FALSE";
    }

    public String getDecoyStatus() {
        return isReverseHit? "TRUE" : "FALSE";
        //return isReverseHit? isDecoyHit() : "FALSE";
        //return isDecoyHit();
    }

    public String getSpectrumCountByExperiment() {
        Peptide [] pipts = getPeptides();
        StringBuffer sb = new StringBuffer();
        for(int i = 0; i < pipts.length; i++) {
            Peptide p = pipts[i];
            int speccount = 0;
            if(p != null) {
                speccount = p.getRedundancyValue();
            } 
            sb.append(speccount + "\t");
        } 
        return sb.toString();
    }

    public Peptide getBestPeptideHit() {
        Peptide bestpep = peptides[0];
        for(int i = 1; i < peptides.length; i++) {
            Peptide p = peptides[i];
            if(p != null) {
                if(bestpep == null) {
                    bestpep = p;
                } else {
                    if(p.getXCorrValue() > bestpep.getXCorrValue()) {
                        bestpep = p;
                    }                   
                }
            }
        }

        return bestpep;

    }

    public PeptideItem(String sequence, int chargestate, boolean isreversehit, int scoretype, int numexp) {
        seq = sequence;
        charge = chargestate;
        isReverseHit = isreversehit;
        scoreType = scoretype;
        peptides = new Peptide[numexp]; 
//System.out.println("Number of experiments: " + numexp);
    }

    public int getNumProteinGroups() {
        return protGroups.size();
    }

    public Iterator<ProteinGroup> getNumProteinGroupMembers() {
        return protGroups.iterator();
    }
    public void assignProteinGroup(ProteinGroup pg) {
        protGroups.add(pg);
    }

    public Peptide [] getPeptides() {

        return peptides;
    }
    public int getNumExperiment() {
        return peptides.length;
    }
    public void addConfidenceScore(double d) {
        confidenceScores.add(d);
    }

    public double getConfidenceProduct() {
        //return xcorrConfidence * zscoreConfidence;
        double product = 1;
        for(int i = 0 ; i < confidenceScores.size(); i++) {
            product *= confidenceScores.get(i);
        }
        return product;
    }

    public double getConfidenceSum() {
        //return xcorrConfidence + zscoreConfidence;
        double sum = 1;
        for(int i = 0 ; i < confidenceScores.size(); i++) {
            sum += confidenceScores.get(i);
        }
        return sum;
    }
    public void setXcorrConfidence(double d) {
        xcorrConfidence = d;
    }
    public double getXcorrConfidence() {
        return xcorrConfidence;
    }

    public void setZscoreConfidence(double d) {
        zscoreConfidence = d;
    }
    public double getZscoreConfidence() {
        return zscoreConfidence;
    }
    public void setConfidence(double d) {
        confidence = d;
    }
    public double getConfidence() {
        return confidence;
    }
    public int compareTo(Object o) {
        PeptideItem pi = (PeptideItem)o;
//        double num1 = getSumZScore();
//        double num2 = pi.getSumZScore(); 

        double num1 = getBestZScore();
        double num2 = pi.getBestZScore(); 
        switch(scoreType) {
            
            case 1 : num1 = getBestZScore(); num2 = pi.getBestZScore(); break;
            case 2 : num1 = getSumZScore(); num2 = pi.getSumZScore(); break;
            case 3 : num1 = getSumXCorr(); num2 = pi.getSumXCorr(); break;
            case 4 : num1 = getConfidence(); num2 = pi.getConfidence(); break;
            case 6 : num1 = getConfidenceSum(); num2 = pi.getConfidenceSum(); break;
            case 7 : num1 = getConfidenceProduct(); num2 = pi.getConfidenceProduct(); break;
            default : num1 = getSumZScore(); num2 = pi.getSumZScore(); break;
        }

        if(num1 == num2) {
            return 0;
        } else if(num1 < num2) {
            return 1;
        } else {
            return -1;
        }
    }
    public boolean containsBothForwardAndReversedHits() {
        if(proteins.size() == 0) {
            return false;
        } 
        Iterator<ProteinItem> it = proteins.iterator();
        ProteinItem pi = it.next();
        boolean isReversedHit = pi.getFasta().isReversed();
        while(it.hasNext()) {
            pi = it.next();
            if(pi.getFasta().isReversed() != isReversedHit) {
                return true;
            } 
        }
        return false;
    }
    public String getSortingCode() {
        StringBuffer sb = new StringBuffer();
        for(int i = 0; i < peptides.length; i++) {
            if(peptides[i] == null) {
                sb.append("0_");
            } else {
                sb.append("1_");
            }
        }
        return sb.toString();
    }
    public void isReverseHit(boolean ir) {
        isReverseHit = ir;
    }
    public boolean isReverseHit() {
        return isReverseHit;
    }   
    public void setFdr(double f) {
        fdr = f;
    } 
    public double getFdr() {
        return fdr;
    } 
    public void addPeptide(Peptide p, int exp) {
        //peptides.add(exp, p);
        peptides[exp] = p;
    }

    public Peptide getPeptide(int expindex) {
        //peptides.get(expindex);
        return peptides[expindex];
    }
    public void addProteinItem(ProteinItem p) {
        proteins.add(p);
    }
    public double getCompositeScore() {
        double score = 0;

        return score;

    }
    public int getCharge() {
        return charge;
    }
    public String getSequence() {
        return seq;
    }
    public String getSeqWithoutMod() {
        //return peptides.get(0).getSeqWithNoModification();
        for(int i = 0; i < peptides.length; i++) {
            if(peptides[i] != null) {

                return peptides[i].getSeqWithNoModification();
            }
        }
        return null;
    }
    public double getBestZScore() {
        if(bestZScore == 0) getSumScores(); 
  
        return bestZScore;
    }  
    public double getBestXCorr() {
        if(bestXCorr == 0) getSumScores(); 
  
        return bestXCorr;
    }  
    public double getSumZScore() {
        if(sumZScore == 0) getSumScores(); 
  
        return sumZScore;
    }  

    public double getAvgZScore() {
        if(sumZScore == 0) getSumScores(); 
  
        return sumZScore/getOccurrence();
    }  
    public double getScore(int scoretype) {
        switch(scoretype) {
            case 1 : return getBestZScore(); 
            case 2 : return getSumZScore(); 
            case 3 : return getSumXCorr(); 
            default : return getSumZScore();
        } 
    }  
 
 
    public double getScore() {
        switch(scoreType) {
            case 1 : return getBestZScore(); 
            case 2 : return getSumZScore(); 
            case 3 : return getSumXCorr(); 
            default : return getSumZScore();
        } 
    }  

    public double getSumXCorr() {
        if(sumXCorr == 0) getSumScores(); 
  
        return sumXCorr;
    }  
    public double getDeltaCn() {
        if(sumDeltaCn == 0) getSumScores(); 
  
        return sumDeltaCn;
    }  
    public int getOccurrence() {
        if(occurance != 0) return occurance;

        for(int i = 0; i < peptides.length; i++) {
           Peptide p = peptides[i];
           if(p != null) {
               occurance++;
           }
         
        }
        return occurance;
    }  
    public int getTotalSpectrumCount() {
        if(totalSpectrumCount == 0) {
            for(int i = 0; i < peptides.length; i++) {
                Peptide p = peptides[i];
                if(p != null) {
                    totalSpectrumCount += p.getRedundancyValue(); 
                }
            
            }
        }
        return totalSpectrumCount;
    }  
    
    public void increaseSpectrumCount() {
        
            for(int i = 0; i < peptides.length; i++) {
                Peptide p = peptides[i];
                if(p != null) {
                     p.getRedundancyValue(); 
                }
            
            }
    }  
    
    public void setTotalSpectrumCount(int value)
    {
        totalSpectrumCount = value;
    }
    
    
    
    private void getSumScores() {
        
        for(int i = 0; i < peptides.length; i++) {
           Peptide p = peptides[i];
           if(p != null) {
                int speccount = p.getRedundancyValue(); 
                //sumDeltaCn += p.getDeltCnValue()*speccount;
                //sumXCorr += p.getXCorrValue()*speccount;
                //sumZScore += p.getSpScoreValue()*speccount;           
               
                double xcorr = p.getXCorrValue();
                sumDeltaCn += p.getDeltCnValue();
                sumXCorr += xcorr;

                double zscore = p.getSpScoreValue();
                sumZScore += zscore;
                if(zscore > bestZScore) {
                    bestZScore = zscore;
                }        
              
                if(xcorr > bestXCorr) {
                    bestXCorr = xcorr;
                }        
                
           }
         
        }
    }  
}

