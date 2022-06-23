/**
 * @file ProteinDatabase.java
 * This is the source file for edu.scripps.pms.util.spectrum.ProteinDatabase
 * @author Tao Xu
 * @date $Date: 2013/01/08 00:36:00 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;
import java.io.FileInputStream;
import java.io.File;
import java.io.IOException;
import java.io.*;

public class ProteinDatabase implements Searchable {

    protected int maxMisCleavage = -1; // -1 means unlimited miscleavage allowed
    protected Protease protease;
    protected int enzymeSpecificity = 0;
    protected ArrayList<Fasta> sequences = new ArrayList<Fasta>(100000);
    //protected HashMap<String, Fasta> ac2Fasta = new HashMap<String, Fasta>(100000);
    protected HashMap<String, Fasta> ac2Fasta = null;
    
    private int[][] freq = new int[5][41];   
    //private boolean isDeCharged = false; 
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
            //System.out.println("ac: " + f.getAccession() + "\tdefline: " + f.getDefline());
        }
      //  //System.out.println("Number of proteins in the database: " + i);
        fis.close(); 
        //fastas = null; // let gc remove the Iterator object
    }   
    public void setMaxMisCleavage(int mc) {
        maxMisCleavage = mc;
    }
    public void setProtease(Protease p) {
        protease = p;
    }
    public void setEnzymeSpecificity(int es) {
        enzymeSpecificity = es;
    }
    protected void addPeptideHit(SearchResult sr, Fasta f, int start, int end) {
        if(enzymeSpecificity == 0 || (protease.checkEnzymeSpecificityStrict(f, start, end) >= enzymeSpecificity)) {
            if(maxMisCleavage == -1 || (protease.getNumInternalMissCleavage(f, start, end) <= maxMisCleavage)) {
                sr.addPeptideHit(f, start, end);
            }
        }
    }
    public static boolean checkTermDiffMods(Fasta f, int start, int end, Modifications m)
    {
        boolean hasCorrectCTerm = true;
        boolean hasCorrectNTerm = true;
        String fastaSeq = f.getSequence();
        char peptideNTerm = fastaSeq.charAt(start);
        char peptideCTerm = fastaSeq.charAt(end);

        if(m.getNTerm()!=null  && m.getNTerm().getSymbol()!= '*')
        {
            hasCorrectNTerm = m.getNTerm().getSymbol() == peptideNTerm;
        }
        if(m.getCTerm() !=null  && m.getCTerm().getSymbol()!= '*' )
        {
            hasCorrectCTerm = m.getCTerm().getSymbol() == peptideCTerm;
        }
        return hasCorrectCTerm && hasCorrectNTerm;
    }

    protected void addPeptideHit(SearchResult sr, Fasta f, int start, int end, Modifications m) {

        if(enzymeSpecificity == 0 || (protease.checkEnzymeSpecificityStrict(f, start, end) >= enzymeSpecificity)) {
            if(maxMisCleavage == -1 || (protease.getNumInternalMissCleavage(f, start, end) <= maxMisCleavage)) {

                if(checkTermDiffMods(f,start,end,m))
                {
                    sr.addPeptideHit(f, start, end, m);
                }

            }
        }
    }

    private void populateAc2Fasta() {
        ac2Fasta = new HashMap<String, Fasta>(1000000);
        for(Iterator<Fasta> it = sequences.iterator(); it.hasNext();) {
            Fasta f = it.next();
            ac2Fasta.put(f.getAccession(), f);
        }
    }
    public Fasta accession2Fasta(String ac) {
        if(ac2Fasta == null) {
            populateAc2Fasta();
        }
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
   
    /* 
    protected PeptideHit createPeptideHit(Peptide p) {
        if(isDeCharged) {
            return new DeChargedPeptideHit(p);
        } else {
           return new PeptideHit(p);
        }
    }
    protected ModifiedPeptideHit createModifiedPeptideHit(Peptide p, Modifications m) {
        if(isDeCharged) {
            return new DeChargedModifiedPeptideHit(p, m);
        } else {
           return new ModifiedPeptideHit(p, m);
        }
    }
    */
    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit)  {
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
                    //result.addPeptideHit(createPeptideHit(new Peptide(f, templeft, rightEnd)));
                    //result.addPeptideHit(f, templeft, rightEnd);
                    addPeptideHit(result, f, templeft, rightEnd);
                } 

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
    
                    addPeptideHit(result, f, leftEnd, rightEnd);
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
                //if (tempmass <= highLimit && rightEnd - templeft > params.getMinimumPeptideLength()-2) {
                if (tempmass <= highLimit) {  
                    //result.addPeptideHit(new ModifiedPeptideHit(f, templeft, rightEnd, m));
                    //result.addPeptideHit(createModifiedPeptideHit(new Peptide(f, templeft, rightEnd), m));
                    //result.addPeptideHit(f, templeft, rightEnd, m);
                    addPeptideHit(result, f, templeft, rightEnd, m);
                }

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
    
                    addPeptideHit(result, f, leftEnd, rightEnd, m);
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
                            //result.addPeptideHit(new ModifiedPeptideHit(f, templeft, rightEnd, m));
                            //result.addPeptideHit(createModifiedPeptideHit(new Peptide(f, templeft, rightEnd), m));
                            //result.addPeptideHit(f, templeft, rightEnd, m);
                            addPeptideHit(result, f, templeft, rightEnd, m);
                            break;
                        }
                    }
                } 

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
                    for(int i = 0; i < numPeaks; i++) {
                        //if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
                        if(mass >= lowLimits[i] && mass <= highLimits[i]) {
                            addPeptideHit(result, f, leftEnd, rightEnd, m);
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
            
        if(numIsotopes == 0) { // for low resolution data, traditional sequest like 
            double highLimit = Math.round(prcMass + ppl.getSearchParams().getHighPrecursorTolerance()/1000);
            //double lowLimit = Math.round(prcMass - ppl.getSearchParams().getHighPrecursorTolerance()/1000);
            double lowLimit = Math.round(prcMass - ppl.getSearchParams().getLowPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit, m);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data
            double diffs = prcMass*acc;
            double highLimit = prcMass + diffs;
            double lowLimit = prcMass - diffs;
            getPeptideHits(sr, highLimit, lowLimit, m);
        } else { // numIsotopes >= 2, for non deisotoped high resolution data

            double [] highLimits = new double[numIsotopes]; 
            double [] lowLimits = new double[numIsotopes]; 
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                highLimits[i] = lowLimits[i] + diffs; 
            }
            getPeptideHits(sr, highLimits, lowLimits, numIsotopes, m);
        } 
           
    }
//    public void getPeptideHits(SearchResult sr, MassCalculator mc, double highLimit, double lowLimit) {
    // precMassAccuracy in ppm
    public SearchResult search(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        precMasses = ppl.getMassCalculator().getPrecMasses();
        //isDeCharged = ppl.isDeCharged();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();
        //MassCalculator mc = ppl.getMassCalculator();
 
        //double mass = prcMass; 
          
        if(numIsotopes == 0) { // for low resolution, traditional sequest like  
            double highLimit = (prcMass + params.getHighPrecursorTolerance()/1000);
            double lowLimit = (prcMass - params.getLowPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data 
            double diffs = prcMass*acc;
            double highLimit = prcMass + diffs;
            double lowLimit = prcMass - diffs;
            getPeptideHits(sr, highLimit, lowLimit);
        } else { //for non deisotoped high resolution data
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
            //System.out.println("finished diffmodSearch\n\n\Zn\n");
        }
        */
        double minPrecMass = ppl.getSearchParams().getMinPrecursorMass();
        double maxPrecMass = ppl.getSearchParams().getMaxPrecursorMass();
        
        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next(); 
            if(m != null && m.getDiffModsShift() != 0 ) {
 //System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                double prcmass = ppl.getPrecursorMass() - m.getMassShift();
                if(prcmass >= minPrecMass && prcmass <= maxPrecMass) {
                    diffModSearch(ppl, sr, m);
                }
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
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
//if(mass > 2040.070 && mass < 2040.072)
//System.out.println("Found!!!! mass: " + mass);
                            //result.addPeptideHit(new PeptideHit(f, templeft, rightEnd));
                            //result.addPeptideHit(createPeptideHit(new Peptide(f, templeft, rightEnd)));
                            //result.addPeptideHit(f, templeft, rightEnd);
                            addPeptideHit(result, f, templeft, rightEnd);
                            break;
                        }
                    }
                }

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
                    for(int i = 0; i < numPeaks; i++) {
                        //if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
                        if(mass >= lowLimits[i] && mass <= highLimits[i]) {
                            addPeptideHit(result, f, leftEnd, rightEnd);
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
    protected void countIsolucineProteins() {
        int numI = 0;
        int numNoI = 0;
        for(Fasta f: sequences) {
            if(f.getSequence().indexOf("I") == -1) {
                numNoI++;
                System.out.println(f.getLength());
            } else {
                numI++;
            }
        }

        System.out.println("NumWithI " + numI + "\tNumNoI: " + numNoI);
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
         //pd.getTrypticPeptideDristribution();
         pd.countIsolucineProteins();

        for(Iterator<Fasta> it = pd.getFastas(); it.hasNext();) {
            Fasta f = it.next();
            String ac = f.getAccession();
            System.out.println("Seq from f: \n" + f.getSequence());
          //  System.out.println("Seq from hash: \n" + pd.accession2Fasta(ac).getSequence());
        } 

    }
}
