package edu.scripps.pms.util.dtaselect;

/**
 * @author Tao Xu 
 * @version $Id
 */
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

public class ProteinItem
{
    Protein protein;
    private HashSet<String> peptides = new HashSet<String>(1000); 
    private int freq = 1; 
    private boolean isProblematic = false;
    private double seqCov = 0;
    private int numTypticTermini = 0;
    private double totalPeptideConfidence = 0;
    private int totalPeptideFrequence = 0;

    public ProteinItem(Protein p) {
        protein = p;
        addProtein(p);
    }
    public double getCompositeScore() {
//System.out.println("seqCov: " + getSeqCoverage()); 
        double confi = getAvgPeptideConfidence()/100;
        confi *= confi;
        confi *= confi;
        confi *= confi;
        confi *= confi;
        return freq*getSeqCoverage()*getNumUniquePeptides() *confi * getNumTrypticTermini() * getAvgPeptideFrequence();///getLength();
    }
    public static String getHeaderLine() {
       
        //return "Accession\tCompositeScore\tFrequence\tnumUniquePeptides\tpeptideFreq\tavgPeptideFreq\t" +
        //       "peptideConf\tavgPeptideConf\tSequenceCoverage\tAvgTypticTermin\tNumTrypticTermini\tProteinLength";
        return "Accession\tFrequence\tnumUniquePeptides\tavgPeptideFreq\t" +
               "avgPeptideConf\tSequenceCoverage\tAvgTypticTermin\tProteinLength";
    }
    public String output() {
       /*return getAccession() + "\t" + getCompositeScore() + "\t" + getFrequence() +
                     "\t" + getNumUniquePeptides() + "\t" + getTotalPeptideFrequence() +
                     "\t" + getAvgPeptideFrequence() + "\t" + getTotalPeptideConfidence() + 
                     "\t" + getAvgPeptideConfidence() + "\t" + getSeqCoverage() +
                     "\t" +  getAvgTrypticTermini() + "\t" + getNumTrypticTermini() + "\t" + getLength();
      */
       return getAccession() + "\t" + getFrequence() +
                     "\t" + getNumUniquePeptides() +  
                     "\t" + getAvgPeptideFrequence() +  
                     "\t" + getAvgPeptideConfidence() + "\t" + getSeqCoverage() +
                     "\t" +  getAvgTrypticTermini() + "\t" + getLength();
    }
    public double getTotalPeptideFrequence() {
        return totalPeptideFrequence;
    }
    public double getAvgPeptideFrequence() {
        return (totalPeptideFrequence + 0.0)/getNumUniquePeptides();
    }
    public double getAvgPeptideConfidence() {
        return totalPeptideConfidence/totalPeptideFrequence;
    }
    public double getTotalPeptideConfidence() {
        return totalPeptideConfidence;
    }
    public String getAccession() {
        return  protein.getLocus();
    }
    public double getAvgTrypticTermini() {
        return getNumTrypticTermini()/(getTotalPeptideFrequence());
    }
    public int getNumTrypticTermini() {
        if(numTypticTermini == 0) {
            for(String peptide : peptides) {
                numTypticTermini += numTrypticTermini(peptide);
            }
        }
        return numTypticTermini;
    }
    
    private static int numTrypticTermini(String peptide) {
        int num = 0;
        char nterm = peptide.charAt(0);
        char cterm = peptide.charAt(peptide.length()-3);
        if(nterm == '-' || nterm == 'K' || nterm == 'R') {
            num++;
        }
        if(cterm == '-' || cterm == 'K' || cterm == 'R') {
            num++;
        }
        return num;
    }
    public void addProtein(Protein p) {
        for(Iterator<Peptide> it = p.getPeptides(); it.hasNext();) {
            Peptide pep = it.next(); 
            String seq = pep.getSequence();
            numTypticTermini  += numTrypticTermini(seq); 
            //peptides.add(seq + pep.getChargeState());
            peptides.add(seq);
            double conf =  Double.parseDouble(pep.getConf());
            //if(conf < 0.8) {
            //    System.out.println("Confidence: " + pep.getConf() + "\t" + conf + "\t" + seq);
            //}
            totalPeptideConfidence += conf;
            totalPeptideFrequence++;
        }
        isProblematic |= p.isProblematic();
        double cov = 0;
        try {
            cov = Double.parseDouble(p.getSeqCoverage().split("%")[0]);
        } catch (Exception e) {
        }
        seqCov = seqCov > cov? seqCov : cov;
//System.out.println("cov: " + cov + "seqCov: " + seqCov);
        freq++;
    }
    public boolean isProblematic() {
        return isProblematic;
    }
    public int getNumUniquePeptides() {
        return peptides.size();
    }
    public int getFrequence() {
        return freq;
    }
    public Iterator<Peptide> getPeptides()
    {
        return  protein.getPeptides();
    }

     

    public String getLocus()
    {
        return  protein.getLocus();
    }

    public String getSeqCount()
    {
        return  protein.getSeqCount();
    }

    public String getSpectrumCount()
    {
        return  protein.getSpectrumCount();
    }

    public double getSeqCoverage()
    {
        return seqCov; 
    }

    public int getLength()
    {
        return  Integer.parseInt(protein.getLength());
    }

    public String getMolWt()
    {
        return  protein.getMolWt();
    }

    public String getPI()
    {
        return  protein.getPI();
    }

    public String getValidation()
    {
        return  protein.getValidation();
    }

    public String getDescription() {
        return  protein.getDescription();
    }

    public int getPeptideSize()
    {
        return  protein.getPeptideSize();
    }


}
