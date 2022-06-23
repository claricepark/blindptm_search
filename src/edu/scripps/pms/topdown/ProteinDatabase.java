/**
 * @file ProteinDatabase.java
 * This is the source file for edu.scripps.pms.util.spectrum.ProteinDatabase
 * @author Tao Xu
 * @date $Date: 2009/08/01 00:51:55 $
 */



package edu.scripps.pms.topdown;

//import edu.scripps.pms.util.seq.TopdownProtein;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;
import java.io.FileInputStream;
import java.io.File;
import java.io.IOException;
import java.io.*;

public class ProteinDatabase {

    protected ArrayList<TopdownProtein> sequences = new ArrayList<TopdownProtein>(100000);
    protected HashMap<String, TopdownProtein> ac2TopdownProtein = new HashMap<String, TopdownProtein>(100000);
    protected SearchParams params; 
    private double [] precMasses;
    // databaseName - the path and name of the protein database
    public ProteinDatabase(String databaseName, MassCalculator mc, SearchParams sp) throws IOException {
        params = sp;
        FileInputStream fis = new FileInputStream(new File(databaseName));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis); 
        int i = 0;
        int maxProteinLength = sp.getMaximumProteinLength();
        System.out.println("Protein database: " + databaseName);
        System.out.println("Maximum protein length: " + maxProteinLength);
        int numseq = 0;
        while(fastas.hasNext()) {
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            numseq++;
            if(f.getLength() <= maxProteinLength) { 
                TopdownProtein tp = new TopdownProtein(f, sp, mc);
                sequences.add(i++, tp);
                ac2TopdownProtein.put(tp.getAccession(), tp);
                System.out.print(i + "\r");
            }
        }
        System.out.println("\nNumber of proteins in the fasta database: " + numseq);
        System.out.println("Number of proteins to be searched: " + i);
        fis.close(); 
        //fastas = null; // let gc remove the Iterator object
    }   
    public TopdownProtein accession2TopdownProtein(String ac) {
        return ac2TopdownProtein.get(ac);        
    }
    public int getNumSequences() {
        return sequences.size();
    }
    public Iterator<TopdownProtein> getTopdownProteins() {
        return sequences.iterator();
    }
    public TopdownProtein getTopdownProtein(int index) {
        return sequences.get(index);
    }

    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit) {
        for(TopdownProtein f : sequences) {
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    break;
                }
                // now the mass should be >= lowLimit
                if (mass <= highLimit) {  
                    result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                    mass -= precMasses[seq[leftEnd++]]; 
                } else {
                    while (mass > highLimit && leftEnd < rightEnd) {
                        mass -= precMasses[seq[leftEnd++]]; 
                    }
                    if (mass >= lowLimit) {
                        result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                        mass -= precMasses[seq[leftEnd++]];
                    }
                }
            }
        }
    }
    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit, Modifications m) {
        for(TopdownProtein f : sequences) {
//System.out.println("in diff");
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    break;
                }
                // now the mass should be >= lowLimit
                if (mass <= highLimit) {  
                    result.addPeptideHit(new ModifiedPeptideHit(f, leftEnd, rightEnd, m));
                    mass -= precMasses[seq[leftEnd++]]; 
                } else {
                    while (mass > highLimit && leftEnd < rightEnd) {
                        mass -= precMasses[seq[leftEnd++]]; 
                    }
                    if (mass >= lowLimit) {
                        result.addPeptideHit(new ModifiedPeptideHit(f, leftEnd, rightEnd, m));
                        mass -= precMasses[seq[leftEnd++]];
                    }
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
        for(TopdownProtein f : sequences) {
           //num++;
        //System.out.println("numproteins: " + num);
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
         // mass of H2O need to be added 
         // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
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
                // now the mass should be >= lowLimit
                if (mass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(mass >= lowLimits[i] && mass <= highLimits[i]) {
                            result.addPeptideHit(new ModifiedPeptideHit(f, leftEnd, rightEnd, m));
                            break;
                        }
                    }
                    mass -= precMasses[seq[leftEnd++]]; 
                } else {
                    while (mass > highLimit && leftEnd < rightEnd) {
                        mass -= precMasses[seq[leftEnd++]]; 
                    }
                    if (mass >= lowLimit) {
                        for(int i = 0; i < numPeaks; i++) {
                            if(mass >= lowLimits[i] && mass <= highLimits[i]) {
                                result.addPeptideHit(new ModifiedPeptideHit(f, leftEnd, rightEnd, m));
                                break;
                            }
                        }
                        mass -= precMasses[seq[leftEnd++]];
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
System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
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
        for(TopdownProtein f : sequences) {
           //num++;
        //System.out.println("numproteins: " + num);
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
         // mass of H2O need to be added 
         // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
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
                // now the mass should be >= lowLimit
                if (mass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(mass >= lowLimits[i] && mass <= highLimits[i]) {
//if(mass > 2040.070 && mass < 2040.072)
//System.out.println("Found!!!! mass: " + mass);
                            result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                            break;
                        }
                    }
                    mass -= precMasses[seq[leftEnd++]]; 
                } else {
                    while (mass > highLimit && leftEnd < rightEnd) {
                        mass -= precMasses[seq[leftEnd++]]; 
                    }
                    if (mass >= lowLimit) {
                        for(int i = 0; i < numPeaks; i++) {
                            if(mass >= lowLimits[i] && mass <= highLimits[i]) {
//if(mass > 2040.070 && mass < 2040.072)
//System.out.println("Found!!!! mass: " + mass);
                                result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                                break;
                            }
                        }
                        mass -= precMasses[seq[leftEnd++]];
                    }
                }
            }
        }
    }
    public SearchResult topDownSearch(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        for(TopdownProtein f : sequences) {
            sr.addPeptideHit(new PeptideHit(f, 0, f.getLength() - 1));
        }
        //MassCalculator mc = ppl.getMassCalculator();
 
        //double mass = prcMass; 
          
        return sr;
    }
    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        SearchParams sp = new SearchParams("topdown.xml");
        ProteinDatabase pd = new ProteinDatabase(sp.getDbName(), new MassCalculator(sp), sp);       
        for(Iterator<TopdownProtein> it = pd.getTopdownProteins(); it.hasNext();) {
            TopdownProtein f = it.next();
            String ac = f.getAccession();
            //System.out.println("Seq from f: \n" + f.getSequence());
            //System.out.println("Seq from hash: \n" + pd.accession2TopdownProtein(ac).getSequence());
        }
          
    }
}
