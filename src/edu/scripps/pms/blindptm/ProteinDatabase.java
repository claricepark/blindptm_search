/**
 * @file ProteinDatabase.java
 * This is the source file for edu.scripps.pms.blindptm.ProteinDatabase
 * @author Tao Xu
 * @date $Date: 2007/09/24 21:30:27 $
 */



package edu.scripps.pms.blindptm;

import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;
import java.io.FileInputStream;
import java.io.File;
import java.io.IOException;
import java.io.*;

public class ProteinDatabase implements Searchable {

    protected ArrayList<Fasta> sequences = new ArrayList<Fasta>(100000);
    protected HashMap<String, Fasta> ac2Fasta = new HashMap<String, Fasta>(100000);
    
    private int[][] freq = new int[5][41];   
    private boolean isDeCharged = false; 
    private double [] precMasses;
    // databaseName - the path and name of the protein database
    public ProteinDatabase(String databaseName) throws IOException {
        FileInputStream fis = new FileInputStream(new File(databaseName));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis); 
        int i = 0;
        while(fastas.hasNext()) {
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            sequences.add(i++, f);
            ac2Fasta.put(f.getAccession(), f);
        }
      //  //System.out.println("Number of proteins in the database: " + i);
        fis.close(); 
        //fastas = null; // let gc remove the Iterator object
    }   
    public Fasta accession2Fasta(String ac) {
        return ac2Fasta.get(ac);        
    }
    public int getNumSequences() {
        return sequences.size();
    }
    public Iterator<Fasta> getFastas() {
        return sequences.iterator();
    }
    public Fasta getFasta(int index) {
        return sequences.get(index);
    }
    
    protected PeptideHit createPeptideHit(Peptide p) {
        if(isDeCharged) {
            return new DeChargedPeptideHit(p);
        } else {
           return new PeptideHit(p);
        }
    }
    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit) {
        for(Fasta f : sequences) {
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0; 
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                while(mass == tempmass) { // ignore non-AA characters
                    mass -= precMasses[seq[leftEnd++]]; 
                }
               
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    //result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                    result.addPeptideHit(createPeptideHit(new Peptide(f, templeft, rightEnd)));
                } 
            }
        }
    }
    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit, Modifications m) {
        for(Fasta f : sequences) {
//System.out.println("in diff");
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0;
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    result.addPeptideHit(new ModifiedPeptideHit(f, templeft, rightEnd, m));
                }
            }
        }
    }
    protected void getPeptideHits(SearchResult result, double [] highLimits, 
                        double [] lowLimits, int numPeaks, Modifications m) {
        int num = 0; 
        double lowLimit = lowLimits[numPeaks-1];
        double highLimit = highLimits[0];
//System.out.println("low: " + lowLimit + "\thigh: " + highLimit);
        for(Fasta f : sequences) {
           //num++;
        //System.out.println("numproteins: " + num);
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
         // mass of H2O need to be added 
         // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0; 
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    //return;
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }

                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
                            result.addPeptideHit(new ModifiedPeptideHit(f, templeft, rightEnd, m));
                            break;
                        }
                    }
                } 
            }
        }
    }
    private void diffModSearch(ProcessedPeakList ppl, SearchResult sr, Modifications m) throws IOException {
        //Modifications m = ppl.getModifications();
        // temp solution, will use m.getMassShift()
        //double massShift = m.getDiffNTermMod();
        double massShift = m.getMassShift();

        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - massShift;
        //MassCalculator mc = ppl.getMassCalculator();
            
        if(numIsotopes == 0) { // traditional sequest like 
            double highLimit = Math.round(prcMass + ppl.getSearchParams().getHighPrecursorTolerance()/1000);
            double lowLimit = Math.round(prcMass - ppl.getSearchParams().getHighPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit, m);
        } else {

            double [] highLimits = new double[numIsotopes]; 
            double [] lowLimits = new double[numIsotopes]; 
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                highLimits[i] = lowLimits[i] + diffs; 
            }
            getPeptideHits(sr, highLimits, lowLimits, numIsotopes, m);
        } 
//System.out.println("HighLimit: " + highLimit + "\tlowLimit: " + lowLimit);
//          System.out.println("diff: " + diff + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);
           
    }
//    public void getPeptideHits(SearchResult sr, MassCalculator mc, double highLimit, double lowLimit) {
    // precMassAccuracy in ppm
    public SearchResult search(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        precMasses = ppl.getMassCalculator().getPrecMasses();
        isDeCharged = ppl.isDeCharged();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();
        //MassCalculator mc = ppl.getMassCalculator();
 
        //double mass = prcMass; 
          
        if(numIsotopes == 0) { // traditional sequest like 
            double highLimit = (prcMass + params.getHighPrecursorTolerance()/1000);
            double lowLimit = (prcMass - params.getLowPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit);
        } else {
            double [] highLimits = new double[numIsotopes]; 
            double [] lowLimits = new double[numIsotopes]; 
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                highLimits[i] = lowLimits[i] + diffs; 
            }
            getPeptideHits(sr, highLimits, lowLimits, numIsotopes);
        }
//          System.out.println("diff: " + diff + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);
        /* // old
        Modifications m = ppl.getModifications();
        if(m != null && m.getDiffModsShift() != 0 ) {
            diffModSearch(ppl, sr);
            //System.out.println("finished diffmodSearch\n\n\n\n");
        }
        */
        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next(); 
            if(m != null && m.getDiffModsShift() != 0 ) {
 //System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                diffModSearch(ppl, sr, m);
//System.out.println("finished diffmodSearch\n\n\n\n");
            }
        }
        return sr;
    }

    // working on making the isotopic peak search faster
    protected void getPeptideHits(SearchResult result,  
                             double [] highLimits, double [] lowLimits, int numPeaks) {
        int num = 0; 
        double lowLimit = lowLimits[numPeaks-1];
        double highLimit = highLimits[0];
//System.out.println("low: " + lowLimit + "\thigh: " + highLimit);
//System.out.println("lowLimit: " + lowLimit + "\thighLimit: " + highLimit);
//for(int i = 0; i < numPeaks; i++) {
//System.out.println("lowLimits[" + i + "]: " + lowLimits[i] + "\thighLimits[" + i+ "]: " + highLimits[i]);
//}
        for(Fasta f : sequences) {
           //num++;
        //System.out.println("numproteins: " + num);
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
         // mass of H2O need to be added 
         // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0;
            int rightEnd = -1;
            int leftEnd = 0; 
             
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    //return;
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                if(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
//if(mass > 2040.070 && mass < 2040.072)
//System.out.println("Found!!!! mass: " + mass);
                            result.addPeptideHit(new PeptideHit(f, templeft, rightEnd));
                            break;
                        }
                    }
                }
            }
        }
    }
    public SearchResult topDownSearch(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        for(Fasta f : sequences) {
            sr.addPeptideHit(new PeptideHit(f, 0, f.getLength() - 1));
        }
        //MassCalculator mc = ppl.getMassCalculator();
 
        //double mass = prcMass; 
          
        return sr;
    }
    protected void getTrypticPeptideDristribution() {
        int maxNumResidue = 5;
        int maxLength = 41;
        int[][] freq = new int[maxNumResidue][maxLength];
        int totalI = 0;
        int totalResidue = 0;
        for(Fasta f : sequences) {
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
            totalResidue += seq.length;
            for(byte b : seq) {
                if(b == 'I') {
                    totalI++;
                }
            }
            int rightEnd = 6;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd <= lastIndex) { 
                if(leftEnd == 0 || isTryptic(seq[leftEnd-1])) {
                    rightEnd = leftEnd + 6;
                    while(rightEnd <= lastIndex && rightEnd - leftEnd < maxLength-1) {
                        if(rightEnd == lastIndex || isTryptic(seq[rightEnd])) {
                            int numI = 0;
                            int numMisCleavage = 0;
                            for(int i = leftEnd; i <= rightEnd; i++) {
                                byte b = seq[i];
                                if(b == 'I') {
                                    numI++;
                                }
                                if(b == 'R' || b == 'K') {
                                    numMisCleavage++;
                                }
                                
                            }
                            if(numI > 4) {
                                numI = 4;
                            }
                            if(numMisCleavage < 2) {
                                freq[numI][rightEnd-leftEnd+1]++;
                            }
                        }
                        rightEnd++;
                    }
                    
                }   
                leftEnd++;
            }
        }
        System.out.println("numResidue\t" + totalResidue + "\tnumI\t" + totalI);
        System.out.print("numI/length\t"); 
        for(int i = 7; i < maxLength; i++) {
            System.out.print(i + "\t");
        }
        System.out.println();
        for(int i = 0; i < maxNumResidue; i++) {
            StringBuffer sb = new StringBuffer(500);
            sb.append(i + "\t");
            for(int j = 7; j < maxLength; j++) {
                sb.append(freq[i][j]+ "\t"); 
            }
            System.out.println(sb);
        }
    }
 
    
    private boolean isTryptic(byte b) {
        return b == 'R' || b == 'K';    
    }
    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        ProteinDatabase pd = new ProteinDatabase(args[0]);
         pd.getTrypticPeptideDristribution();
        /*
        for(Iterator<Fasta> it = pd.getFastas(); it.hasNext();) {
            Fasta f = it.next();
            String ac = f.getAccession();
            System.out.println("Seq from f: \n" + f.getSequence());
            System.out.println("Seq from hash: \n" + pd.accession2Fasta(ac).getSequence());
        } 
        */ 
    }
}
