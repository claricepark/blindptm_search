/**
 * @file MemorizedDb.java
 * This is the source file for edu.scripps.pms.mspid.MemorizedDb
 * @author Tao Xu
 * @date $Date: 2005/11/10 00:02:06 $
 */
package edu.scripps.pms.blindptm;

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.ByteArrayConverter;
import java.util.*;
import java.io.*;
//import java.lang.reflect.Array;

// command line processor
import org.apache.commons.cli.*;


class MemorizedDb {
    private static final String USAGE= "---Usage: java MemorizedDb dbfile ---";
    private static final String DBNAMEOPT = "d";
    private static final String TAXONOPT = "t";
    private static final String MOLTYPEOPT = "m";
    private static final String OUTPUTOPT = "o";
  
    /**
     * The protein index start from 0 to number of proteins - 1
     * The .pin the same number of records as the number of peptides,
     * Each record has 4 bytes for mass (equals to Math.round(realMass*ACCURACYFACTOR),
     * 4 bytes for protein index, 2 bytes for start (unsighed short, same as char), 
     * 2 bytes for end (unsighed short, same as char)
     * The .min file contains offset of masses in the .pin file, each record is a long
     * 
     */
    //private TimeUtils timer = new TimeUtils();
    private static long startTime;
    private static long endTime;
    private String fastaFile;
    //private RandomAccessFile raf;
    private static int ACCURACYFACTOR = MassCalculator.MASSACCURACYFACTOR;
    private int LOWLIMIT; 
    private int HIGHLIMIT;
    private final int LENGTHLIMIT = 7;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    private static int NUMBINS;
    // for the frequence of peptide length
    private int [] massFreq;
    private Protease protease; 
    private int enzymeSpecificity; 
    long [] massIndex;
    //private static final int NUMRECORDSPERCHUNK = 10000000;
    private static int maxLength = 0;
    private SearchParams sp = null;
    private MassCalculator mc;
    long totalNumPeptides = 0;
    //private byte [] temp = new byte[RECORDSIZE];
    // the following variable will be used by the heap sort for swapping
    //private String [] fileNames;
    public MemorizedDb(String fasta, SearchParams sp) throws IOException { 
        this.fastaFile = fasta;
        this.sp = sp;
        mc = new MassCalculator(sp);
        protease = sp.getProtease();
        enzymeSpecificity = sp.getEnzymeSpecificity();
        // .pin for peptide index, .min for mass index 

    }
    private void setParameters(String [] args) {

        LOWLIMIT = (int)sp.getMinPrecursorMass()*ACCURACYFACTOR;
        HIGHLIMIT = (int)sp.getMaxPrecursorMass()*ACCURACYFACTOR;
        // keep 1pp accuracy for peptides with mass 1000 dalton
        NUMBINS = HIGHLIMIT + 1; 
        massIndex = new long[NUMBINS];
        // for the frequence of peptide length
        massFreq = new int[NUMBINS];
        // one file for every AMUSPERFILE amu
    }

    private void outputPeptides() throws IOException {
        TimeUtils timer = new TimeUtils();
        timer.startTiming();
        ArrayList<ArrayList<edu.scripps.pms.mspid.Peptide>>  dos = new ArrayList<ArrayList<edu.scripps.pms.mspid.Peptide>>(NUMBINS);
        for(int i = 0; i < NUMBINS; i++) {
            dos.add(new ArrayList<edu.scripps.pms.mspid.Peptide>());
        }

        
        FileInputStream fis = new FileInputStream(new File(fastaFile));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        int i = 0; // index of the fasta
        while(fastas.hasNext()) {
            Fasta f = fastas.next();
            totalNumPeptides += outputPeptides(i++, f, dos);
            System.out.print("Number of proteins processed: " + i + "\r");
        }
        int numZeros = 0;
        int numNonZeros = 0;
        for(int k = LOWLIMIT; k < NUMBINS; k++) {
            if(massFreq[k] > 0) {
                numNonZeros++;
            } else {
                numZeros++;
            }
        }
        timer.stopTiming();
        System.out.println("Number of proteins in the database: " + i + "\tLongest protein: " + maxLength);
        System.out.println("Number of peptides: " + totalNumPeptides + "\tavgNumPeptidesPerProtein: " + totalNumPeptides/(double)i);
        System.out.println("NumZeros: " + numZeros + "\tNumNonZeros: " + numNonZeros);
        System.out.println("Time used to output peptides: " + timer.getTimeUsed());
        System.out.println("Number of proteins in the database: " + i + "\tLongest protein: " + maxLength);
        System.out.println("Number of peptides: " + totalNumPeptides + "\tavgNumPeptidesPerProtein: " + totalNumPeptides/(double)i);
        System.out.println("NumZeros: " + numZeros + "\tNumNonZeros: " + numNonZeros);
        fis.close();
    }

    // index is the index for the fasta in the database
    private int outputPeptides(int index, Fasta f, ArrayList<ArrayList<edu.scripps.pms.mspid.Peptide>>  dos) throws IOException {
        
         // to check the longest protein in the database
        int lastIndex = f.getLength();
        if(lastIndex > maxLength) { 
            maxLength = lastIndex;
        }
        byte [] seq = f.getSequenceAsBytes();
        int numPeptides = 0; 
        for(int i = 0; i < lastIndex; i++) {
            int j = i; // use char as unsigned short
             // mass of H2O need to be added
             // because the mass of the residues does not count H2O
            float mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSH;
            int length = 0;
            while(mass*ACCURACYFACTOR <= HIGHLIMIT && j < lastIndex) {
                
                mass += mc.getPrecursorMass(seq[j]);
                length++;
                int intMass = Math.round(mass * ACCURACYFACTOR);
                if(length >= LENGTHLIMIT && intMass >= LOWLIMIT && intMass <= HIGHLIMIT) {
                    if(protease == null || protease.checkEnzymeSpecificity(f, i, j) >= enzymeSpecificity) {
                        dos.get(intMass).add(new edu.scripps.pms.mspid.Peptide(f, i, j)); // end of the peptide
                        numPeptides++;
                        massFreq[intMass]++; 
                    }
                }
                j++;
            }
        }
        return numPeptides;
    }
    public static void main(String args[]) throws Exception {
        try {
            TimeUtils timer = new TimeUtils();
            SearchParams sp = new SearchParams("search.xml");
            String fastaFile = args[0];
            timer.startTiming();
            MemorizedDb se = new MemorizedDb(fastaFile, sp);
            se.setParameters(args);
            se.outputPeptides();
            // prepare for garbage collection to reduce memory usage
            System.out.println("Finished sorting");
            //se.writeMassIndex();

            timer.stopTiming();
            long timeUsed = timer.getTimeUsedMillis();
            System.out.println("Time used for process database: " + timeUsed);
            System.out.println();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }


}


