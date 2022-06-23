/**
 * @file BlindProteinDatabase.java
 * This is the source file for edu.scripps.pms.blindptm.BlindProteinDatabase
 * @author Tao Xu
 * @date 
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
            int tempright = 0; 
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
//System.out.print("mass: " + mass + "\t" + (char)(seq[rightEnd]) + "\t" + precMasses[seq[rightEnd]] + "\n");
                }
                if (mass < lowLimit) {
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                tempright = rightEnd;
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
               
//System.out.println("templeft: " + templeft + "\ttempright: " + tempright + "\ttempmass: " + tempmass);
                // now the mass should be >= lowLimit
                while (tempmass <= highLimit && tempright < lastIndex) {  


                    double massdiff = prcmass - tempmass;
//System.out.println("templeft: " + templeft + "\ttempright: " + tempright + "\t" + (char)seq[tempright] +"\tempmass: " + tempmass + "\tmassdiff: " + massdiff );
                    Peptide p = new Peptide(f, templeft, tempright);
                    result.addPeptideHit(new PeptideHit(p));
            //result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
//System.out.println("mass: " + tempmass + " massdiff: " + massdiff + "\ttempleft: " + templeft + "\ttempright: " + tempright + " residue: " + (char)seq[tempright]);
                    if(massdiff >= 3 || massdiff <= -0.9) {
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
                    tempmass += precMasses[seq[++tempright]]; 
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
        //double acc = precMassAccuracy/1000000.0f;

        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();

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
        
    }
}
