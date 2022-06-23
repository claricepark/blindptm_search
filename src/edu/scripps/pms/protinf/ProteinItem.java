package edu.scripps.pms.protinf;

import java.util.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.util.dtaselect.Peptide;



// 
public class ProteinItem implements Comparable {
    private HashSet<PeptideItem> peptides = new HashSet();
    private ArrayList<PeptideItem> uniquePeptides = null;
    private double sumZScore = 0;
    private double bestPeptideScore = 0;
    private Fasta fasta = null;
    private int totalSpectrumCount = 0; 
    private double proteinScore = 0;
    private float pi = -1;
    private int numIdentifiedTrypticPeptides = -1;
    private int numTrypticPeptides = -1;
    private ArrayList<Peptide> [] peptidesByExperiment = null;
    private int [] peptideCountsByExperiment = null;
    private int [] spectrumCountsByExperiment = null;

    private double sumUniqueZScore = -1;
    private double bestUniquePeptideScore = -1;
    private double isoformScore = -1;

    private float seqcov = -1;
    private int numuniquepepseq = -1;

    // number of unique peptide sequence ignoring charge and modification status
    public int getNumUniquePeptideSequence() {
        if(numuniquepepseq == -1) {
            HashSet<String> uniqpepseq = new HashSet();
            for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
                PeptideItem pi = it.next();
                String pseq = pi.getSeqWithoutMod();
                uniqpepseq.add(pseq); 
            }

            numuniquepepseq = uniqpepseq.size();
        } 
        return numuniquepepseq;
    }

    public float getSeqCov() {
        if(seqcov == -1) {
            String seq = fasta.getSequence(); 
            boolean [] covered = new boolean[seq.length()];
            for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
                PeptideItem pi = it.next();
                String pseq = pi.getSeqWithoutMod();
                int index = seq.indexOf(pseq); 
                if(index == -1 ) {
System.out.println("Peptide sequence " + pseq + " cannot be found in this protein accession: " + fasta.getAccession() +"\tDescription: " + fasta.getDefline());
                } else {
                    for(int i = 0; i < pseq.length(); i++) {
                        covered[index + i] = true;
                    } 
                }
             
            }
            int post = 0; // postive residue number
            for(int i = 0; i < covered.length; i++) {
                if(covered[i]) {
                    post++;
                }
            }

            if(post == 0 ) {
System.out.println("No peptide cannot be found in this protein " + fasta.getAccession());
            } 
            seqcov = 100*(0.0f + post)/covered.length; 
        } 
        return seqcov;
    }

    public String getSpectrumCountInfo() {
        getSpectrumCountByExperiment();
        StringBuffer sb = new StringBuffer();
        for(int i = 0; i < spectrumCountsByExperiment.length; i++) {
            sb.append(spectrumCountsByExperiment[i]);
            sb.append("\t");
        }
        return sb.toString();
    }
    public String getPeptideCountInfo() {
        getSpectrumCountByExperiment();
        StringBuffer sb = new StringBuffer();
        for(int i = 0; i < peptideCountsByExperiment.length; i++) {
            sb.append(peptideCountsByExperiment[i]);
            sb.append("\t");
        }
        return sb.toString();
    }
    public int [] getSpectrumCountByExperiment() {
        if(spectrumCountsByExperiment == null) {
            calcSpectrumPeptideCountByexperiment();
        }  
         
        return spectrumCountsByExperiment;
    }
    public int [] getPeptideCountByExperiment() {
        if(peptideCountsByExperiment == null) {
            calcSpectrumPeptideCountByexperiment();
        }  
         
        return peptideCountsByExperiment;
    }
    private void calcSpectrumPeptideCountByexperiment() {

        if(peptideCountsByExperiment == null || spectrumCountsByExperiment == null ) {
            getPeptidesByExperiment();
            peptideCountsByExperiment = new int[peptidesByExperiment.length];
            spectrumCountsByExperiment = new int[peptidesByExperiment.length];
            for(int i = 0; i < peptidesByExperiment.length; i++) {
                int speccounts = 0; 
                int peptcounts = 0;
                for(Iterator<Peptide> it = peptidesByExperiment[i].iterator(); it.hasNext();) {
                    Peptide p = it.next();
                    speccounts += p.getRedundancyValue();
                    peptcounts++;

                }

                peptideCountsByExperiment[i] = peptcounts;
                spectrumCountsByExperiment[i] = speccounts;
            }
        }  
         
    }
    public double getSumUniqueZScore() {
        if(sumUniqueZScore != -1 || getUniquePeptideItems().size() < 1) return sumUniqueZScore;
        Iterator<PeptideItem> it =  getUniquePeptideItems().iterator(); 
        sumUniqueZScore = 0;
        while(it.hasNext()) {
            PeptideItem pi = it.next();
            sumUniqueZScore += pi.getSumZScore();
        }
        return sumUniqueZScore;
    }     

    
    public double getBestUniquePeptideScore() {
        if(bestUniquePeptideScore != -1 || getUniquePeptideItems().size() < 1) return bestUniquePeptideScore;
        Iterator<PeptideItem> it =  getUniquePeptideItems().iterator(); 
        while(it.hasNext()) {
            bestUniquePeptideScore += it.next().getScore();
        }
        return bestUniquePeptideScore;
    }     

    public double getIsoformScore() {
        return isoformScore;
    }

    public void setIsoformScore(double ifs) {
        isoformScore = ifs;
    }
    public ArrayList<PeptideItem> getUniquePeptideItems() {
        if(uniquePeptides != null) {
            return uniquePeptides;
        }
        uniquePeptides = new ArrayList();
       for(Iterator<PeptideItem> it = getPeptideItems().iterator(); it.hasNext();) {
           PeptideItem pi = it.next();
           if(pi.getNumProteinGroups() == 1) {
               uniquePeptides.add(pi);
           } else if(pi.getNumProteinGroups() == 0) {
               System.out.println("Something wrong, Found peptideitem with no protein group assgined");
           }
       }
       return uniquePeptides;
    }

    private boolean isNtermTryptic(String protseq, String pepseq) {
        int index = protseq.indexOf(pepseq);
        if(index == 0) {
            return true;
        } else if(index > 0) {

            char c = protseq.charAt(index-1);
            if(c == 'R' || c == 'K') {
                return true;    
            }
        } else {
System.out.println("Weird, pepseq not in protein\t" + pepseq + "\t" + protseq);
        }
        return false;
    }
    public int getNumExperiment() {
        for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
            return it.next().getNumExperiment();
        } 
        return 0;
    }
    public ArrayList<Peptide> [] getPeptidesByExperiment() {
        
        if(peptidesByExperiment == null) {
            peptidesByExperiment = new ArrayList[getNumExperiment()];
            for(int i = 0; i < peptidesByExperiment.length; i++) {
                peptidesByExperiment[i] = new ArrayList();
            }
        } else {
            return peptidesByExperiment;
        }
        for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
            Peptide [] peps = it.next().getPeptides();
            for(int i = 0; i < peps.length; i++) {
                if(peps[i] != null) {
                    peptidesByExperiment[i].add(peps[i]);
                }

            }
        }
        return peptidesByExperiment;
    } 

    private boolean isCtermTryptic(String protseq, String pepseq) {
        int index = protseq.lastIndexOf(pepseq);
        if(index + pepseq.length() == protseq.length()-1) {
            return true;
        } else {

            char c = pepseq.charAt(pepseq.length()-1);
            if(c == 'R' || c == 'K') {
                return true;    
            }
        }
        return false;
    }
    private static int getNumMiscleavage(String pepseq) {
        int nummiscleav = 0;
        for(int i = 0; i < pepseq.length() -2; i++) {
            char c = pepseq.charAt(i);
            if(c == 'K' || c == 'R') {
                nummiscleav++;
            }

        }
        return nummiscleav;
    } 
    public int getNumIdentifiedTrypticPeptides() {
        if(numIdentifiedTrypticPeptides == -1) {
            numIdentifiedTrypticPeptides = 0;
            HashSet<String> uniqseqs = new HashSet();
            String protseq = fasta.getSequence(); 
            for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
                PeptideItem pitem = it.next();
                //uniqseqs.add(pitem.getSequence()); 
                String seqnomod = pitem.getSeqWithoutMod();
                //uniqseqs.add(pitem.getSeqWithoutMod()); 
                uniqseqs.add(seqnomod); 
                int numtrypends = 0;
                if(isNtermTryptic(protseq, seqnomod)) numtrypends++;
                if(isCtermTryptic(protseq, seqnomod)) numtrypends++;
                pitem.setTrypticStatus(numtrypends);  
                
            }
            for(Iterator<String> it = uniqseqs.iterator(); it.hasNext();) {
                String seq = it.next();
                
                if(seq.length() >= 6 && isNtermTryptic(protseq, seq) && isCtermTryptic(protseq, seq) && getNumMiscleavage(seq) == 0) {
                    numIdentifiedTrypticPeptides++;
                } 
            }
        }

        return numIdentifiedTrypticPeptides;
    }

    public int getNumTrypticPeptides(MassCalculator mc) {
        if(numTrypticPeptides == -1) {
            //numTrypticPeptides = ProteinTrypticPeptideNumber.getNumTrypticPeptides(fasta, mc, 4500, 600, 6);
            numTrypticPeptides = ProteinTrypticPeptideNumber.getNumTrypticPeptides(fasta, 6);
 
        }
        if(getNumIdentifiedTrypticPeptides() > numTrypticPeptides) {
            numTrypticPeptides = getNumIdentifiedTrypticPeptides();
        }
        return numTrypticPeptides;
    }
    public int getNumTrypticPeptides() {
        if(numTrypticPeptides == -1) {
            //numTrypticPeptides = ProteinTrypticPeptideNumber.getNumTrypticPeptides(fasta, mc, 4500, 600, 6);
            numTrypticPeptides = ProteinTrypticPeptideNumber.getNumTrypticPeptides(fasta, 6);
 
        }
        if(getNumIdentifiedTrypticPeptides() > numTrypticPeptides) {
            numTrypticPeptides = getNumIdentifiedTrypticPeptides();
        }
        return numTrypticPeptides;
    }
    public float getProteinPI() {
        if(pi == -1) {
            pi = edu.scripps.pms.util.seq.ProteinInfo.calculatePI(getFasta().getSequence())[0];
        }
        return pi;
    }
    public ProteinItem(Fasta f) {
        fasta = f;
    
    }
    public boolean isReverseHit() {
        return fasta.getAccession().startsWith("Revers");
    }
    public Fasta getFasta() {
        return fasta;
    }
    public void addPeptideItem(PeptideItem p) {
        peptides.add(p);
    }
    public HashSet<PeptideItem> getPeptideItems() {
        return peptides;
    }
    public int getTotalSpectrumCount() {
        if(totalSpectrumCount == 0) {
            for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
                PeptideItem p = it.next();
                if(p != null) {
                    totalSpectrumCount += p.getTotalSpectrumCount();
                }

            }
        }
        return totalSpectrumCount;
    }
    public void setProteinScore(double s) {
        proteinScore = s;
    }
    public double getProteinScore() {
        return proteinScore;
    }
    public double getBestPeptideScore() {
       
        if(bestPeptideScore == 0) {
            for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
                //sumZScore += it.next().getSumZScore();
                bestPeptideScore += it.next().getScore();
            }
        }

        return bestPeptideScore;
    }
    public double getSumZScore() {
        if(sumZScore == 0) {
            for(Iterator<PeptideItem> it = peptides.iterator(); it.hasNext();) {
                //sumZScore += it.next().getSumZScore();
                sumZScore += it.next().getScore();
            }
        }

        return sumZScore;
    }

    public double getAverageZScore() {
        return getSumZScore()/getNumPeptides();
    }
    public int getNumPeptides() {
        return peptides.size();
    }

    public int compareTo(Object o) {
        //int num1 = peptides.size();
        //int num2 = ((ProteinItem)o).getPeptideItems().size();
        //double num1 = getSumZScore();
        //double num2 = ((ProteinItem)o).getSumZScore();
        double num1 = getProteinScore();
        double num2 = ((ProteinItem)o).getProteinScore();
        if(num1 == num2) {
            return 0;
        } else if(num1 < num2) {
            return 1;
        } else {
            return -1;
        }

    }

}

