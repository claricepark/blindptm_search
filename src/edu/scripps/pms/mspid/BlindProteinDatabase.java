/**
 * @file BlindProteinDatabase.java
 * This is the source file for edu.scripps.pms.util.spectrum.BlindProteinDatabase
 * @author Tao Xu
 * @date $Date: 2008/01/17 19:46:44 $
 */



package edu.scripps.pms.mspid;

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

public class BlindProteinDatabase implements Searchable {

    protected ArrayList<Fasta> sequences = new ArrayList<Fasta>(100000);
    protected HashMap<String, Fasta> ac2Fasta = new HashMap<String, Fasta>(100000);
   
    private boolean isDeCharged = false; 
    private double [] precMasses;
    private SearchParams parameters;
    private MassCalculator mc;
    // databaseName - the path and name of the protein database
    public BlindProteinDatabase(String databaseName, SearchParams params) throws IOException {
        this.parameters = params;
        mc = new MassCalculator(params);
        FileInputStream fis = new FileInputStream(new File(databaseName));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis); 
        int i = 0;
        while(fastas.hasNext()) {
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            double mplush = mc.getPrecursorMass(f.getSequence());
            f.setMPlusH(mplush);
            sequences.add(i++, f);
            ac2Fasta.put(f.getAccession(), f);
        }
        System.out.println("Number of proteins in the database: " + i);
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
    
    protected void getPeptideHits(SearchResult result, double prcmass) {
        for(Fasta f : sequences) {
//System.out.println("in diff");
            double massdiff = prcmass - f.getMPlusH();
            int lastIndex = f.getLength() - 1;
            Peptide p = new Peptide(f, 0, lastIndex);
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            int rightEnd = -1;
            int leftEnd = 0; 
            result.addPeptideHit(new PeptideHit(p));
            //result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
            if((massdiff >= 6 && massdiff<=100) || (massdiff <= -0.9 && massdiff >= -40)) {

                DiffMod m = new DiffMod(massdiff, '~');
                if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                    result.addPeptideHit(p, m);
                }
                m = new DiffMod(massdiff-MassSpecConstants.MASSDIFFC12C13, '~');
                if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                    result.addPeptideHit(p, m);
                }
                m = new DiffMod(massdiff-2*MassSpecConstants.MASSDIFFC12C13, '~');
                if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                    result.addPeptideHit(p, m);
                }
                m = new DiffMod(massdiff-3*MassSpecConstants.MASSDIFFC12C13, '~');
                if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                    result.addPeptideHit(p, m);
                }
               
            }
        }
    }

    protected void getPeptideHits(SearchResult result, double prcmass, 
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
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
               
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    //result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                    //result.addPeptideHit(createPeptideHit(new Peptide(f, templeft, rightEnd)));

                    double massdiff = prcmass - tempmass;
                    Peptide p = new Peptide(f, templeft, rightEnd);
                    result.addPeptideHit(new PeptideHit(p));
            //result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                    if(massdiff >= 3 || massdiff <= -0.9) {
//System.out.println("massdiff: " + massdiff);
                        DiffMod m = new DiffMod(massdiff, '~');
                        if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                            result.addPeptideHit(p, m);
                        }
                        m = new DiffMod(massdiff-MassSpecConstants.MASSDIFFC12C13, '~');
                        if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                            result.addPeptideHit(p, m);
                        }
                        m = new DiffMod(massdiff-2*MassSpecConstants.MASSDIFFC12C13, '~');
                        if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                            result.addPeptideHit(p, m);
                        }
                        m = new DiffMod(massdiff-3*MassSpecConstants.MASSDIFFC12C13, '~');
                        if(!(m.getMassShift() > -0.9 && m.getMassShift() < 0.9)) {
                            result.addPeptideHit(p, m);
                        }
               
                    }
                } 
            }
        }
    }
    public SearchResult search(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        precMasses = ppl.getMassCalculator().getPrecMasses();
        isDeCharged = ppl.isDeCharged();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;

        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();
        //getPeptideHits(sr, prcMass); // for database with identified peptide only just to prove concept
        double highlimit = prcMass + params.getLowPrecursorTolerance()/1000;
        double lowlimit = prcMass - params.getHighPrecursorTolerance()/1000;
//System.out.println("highlimit: " + highlimit + "\tlowlimit: " + lowlimit);
        getPeptideHits(sr, prcMass, highlimit, lowlimit);
 
          
        return sr;
    }

    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        SearchParams params = new SearchParams("search.xml");
        BlindProteinDatabase pd = new BlindProteinDatabase(params.getDbName(), params);
        double high = 800000;
        double low = 600000;
        int numHigh = 0;
        int numLow = 0;
        for(Iterator<Fasta> it = pd.getFastas(); it.hasNext();) {
            Fasta f = it.next();
            double mass = f.getMPlusH();
            if(mass > high) {
                System.out.println(f.getDefline());
                numHigh++;
            } 
            if(mass < low) {
                numLow++;
            }
            
        }
        
        System.out.println("Number of Proteins in the " + params.getDbName() + " database: " + pd.getNumSequences());
        System.out.println("Number of Proteins with mass >= " + high + ": " + numHigh);
        System.out.println("Number of Proteins with mass < " + low + ": " + numLow);
    }
}
