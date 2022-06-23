package edu.scripps.pms.protinf;

import java.util.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.dtaselect.Peptide;



// 
public class ProteinGroup implements Comparable {

    private ArrayList<ProteinItem> proteins = new ArrayList();
    private ArrayList<ProteinGroup> subset = new ArrayList();
    private ProteinItem representative = null;
    private String accs = null;

    private double sumZScoreConfidence = 0; 
    private double trypticConfidence = 0; // confidence based on identified tryptic peptide fraction
    private double averageZScoreConfidence = 0;
    private double sumConfidence = 0;  
    private double productConfidence = 0;  

    private int proteinGroupId = 0;
    private double finalConfidence = 0;
    private double finalGlobalFdr = 0;
    //private ArrayList<Peptide> [] peptides = null;
    private boolean isAssignedToPeptideItem = false;
    

    public void assignProteinGroupToPeptideItem() {
        if(!isAssignedToPeptideItem) {
            for(Iterator<PeptideItem> it = getRepresentative().getPeptideItems().iterator(); it.hasNext();) {
                PeptideItem pi = it.next();
                pi.assignProteinGroup(this);
            }
            isAssignedToPeptideItem = true;
        }
    }
    public double getIdentifiedTrypticPeptideFraction() {
        if(getRepresentative().getNumTrypticPeptides() <= 0) return 0;
        return getRepresentative().getNumIdentifiedTrypticPeptides()/getRepresentative().getNumTrypticPeptides();
    }

    public String getOccurrenceAndSortingCode() {
        StringBuffer sb = new StringBuffer();
        ArrayList<Peptide> [] peps = getRepresentative().getPeptidesByExperiment(); 
        int occurance = 0;
        for(int i = 0; i < peps.length; i++) {
            if(peps[i].size() > 0) {
                sb.append("1_");
                occurance++;
            } else {
                sb.append("0_");
            }
        }
        return occurance + "\t" + sb.toString();
    }

    public void setProteinGroupId(int id) {
        proteinGroupId = id;
    }
    public void setFinalConfidence(double conf) {
        finalConfidence = conf;
    }
    public void setFinalGlobalFDr(double fdr) {
        finalGlobalFdr = fdr;
    }

    public int getProteinGroupId() {
        return proteinGroupId;
    }
    public double getFinalConfidence() {
        return finalConfidence;
    }
    public double setFinalGlobalFDr() {
        return finalGlobalFdr;
    }

    public void addProteinItem(ProteinItem p) {
        proteins.add(p);
    }
    //public void addSubsetProteinGroup(ProteinGroup pg) {
    //    subset.add(pg);
    //}

    public void setSumZScoreConfidence(double d) {
        sumZScoreConfidence = d;
    }

    public void setTrypticConfidence(double d) {
        trypticConfidence = d;
    }
    public double getTrypticConfidence() {
        return trypticConfidence;
    }
    public void setAverageZScoreConfidence(double d) {
        averageZScoreConfidence = d;
    }

    public double getConfidenceSum() {
        return sumZScoreConfidence + averageZScoreConfidence + 0.2*trypticConfidence;
    }

    public double getConfidenceProduct() {
        return sumZScoreConfidence*averageZScoreConfidence;
    }
    public boolean isReverseHit() {
        for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
            if(!it.next().isReverseHit()) return false;
        }
        return true;
    }
 
    public String getAccs() {
        StringBuffer sb = new StringBuffer(500);    
        for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
            Fasta f = it.next().getFasta();
            sb.append(f.getAccession());
            sb.append(";"); 
        }
        return sb.toString();
    }
    public String getDeflines() {
        StringBuffer sb = new StringBuffer(500);    
        for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
            Fasta f = it.next().getFasta();
            sb.append(f.getDefline());
            sb.append(";"); 
        }
        return sb.toString();
    }
    public int getTotalSpectrumCount() {
        return getRepresentative().getTotalSpectrumCount();
    }
    public String getRepresentativeAcc() {
        return getRepresentative().getFasta().getAccession();
    }
    public int getRepresentativeLength() {
        return getRepresentative().getFasta().getLength();
    }

    public float getRepresentativePI() {
        return getRepresentative().getProteinPI();
    }
    public String getRepresentativeDefline() {
        return getRepresentative().getFasta().getDefline();
    }

    // the shortest protein is assigned as representative protein of this protein group
    public ProteinItem getRepresentative() {
        if(representative == null) {
            for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
                ProteinItem f = it.next();
                if(representative == null) {
                    representative = f;
                } else {
                    if(representative.getFasta().getLength() > f.getFasta().getLength()) {
                        representative = f;
                    }
                } 
            }
        }
        return representative;
    }
    public String getAccessions() {
        StringBuffer sb = new StringBuffer();
        if(accs == null) {
             for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
                 sb.append(it.next().getFasta().getAccession());
                 sb.append(";");
             }
        } 
        return sb.toString();
    }
    public ArrayList<ProteinItem> getProteinItems() {
        return proteins;
    }
    
    public HashSet<PeptideItem> getPeptideItems()  {
        if(proteins.size() > 0) {
            return proteins.get(0).getPeptideItems();
        } else {
            return null;
        }
    }
    public void addSubset(ProteinGroup pg) {
        subset.add(pg);
    }
 
    public int getNumProteinItems() {
        return proteins.size();
    }
    public int getTotalSpectrumCoun() {
        return proteins.get(0).getTotalSpectrumCount();
    }
    public int getNumPeptideItems() {
        return proteins.get(0).getPeptideItems().size();
    }
    public double getSumZScore() {
        return proteins.get(0).getSumZScore();
    }
    public int compareTo(Object o) {
        //int num1 = peptides.size();
        //int num2 = ((ProteinItem)o).getPeptideItems().size();
        //double num1 = getRepresentative().getSumZScore();
        //double num2 = ((ProteinGroup)o).getRepresentative().getSumZScore();
        double num1 = getRepresentative().getProteinScore();
        double num2 = ((ProteinGroup)o).getRepresentative().getProteinScore();
        if(num1 == num2) {
            return 0;
        } else if(num1 < num2) {
            return 1;
        } else {
            return -1;
        }

    }
}

