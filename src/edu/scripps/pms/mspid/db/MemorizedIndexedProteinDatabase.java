/**
 * @file MemorizedIndexedProteinDatabase.java
 * This is the source file for edu.scripps.pms.mspid.MemorizedIndexedProteinDatabase
 * @author Tao Xu
 * @date $Date: 2007/09/24 21:30:27 $
 */
package edu.scripps.pms.mspid.db;

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.ByteArrayConverter;
//import edu.scripps.pms.util.enzyme.Protease;
import java.util.*;
import java.io.*;
import java.lang.reflect.Array;

// assume all peptide index info kept in RAM 
public class MemorizedIndexedProteinDatabase extends ProteinDatabase {
    /**
     * The protein index start from 0 to number of proteins - 1
     * The .pin the same number of records as the number of peptides,
     * Each record has 4 bytes for mass (equals to (int)Math.round(realMass*ACCURACYFACTOR),
     * 4 bytes for protein index, 2 bytes for start (unsighed short, same as char), 
     * 2 bytes for end (unsighed short, same as char)
     * The .min file contains offset of masses in the .pin file, each record is a long
     * 
     */
    private TimeUtils timer = new TimeUtils();
    private static final String USAGE = "java MemorizedIndexedProteinDatabase dbfile";
    private String peptideIndexFile;
    private String massIndexFile;
    //private RandomAccessFile raf;
    private static final int ACCURACYFACTOR = 1000;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    // for the frequence of peptide length
    private static final int RECORDSIZE = 8;
    //private static final int DEFAULTNUMRECORDS = 1000000; 
    private long [] index; // mass index for offsets of each mass
    private int MAXINDEX;

    private int [] proteinIndex;
    private char [] peptideStartIndex;
    private char [] peptideEndIndex;

    long totalNumPeptides = 0;
    // the following variable will be used by the heap sort for swapping
    //private byte [] b = null; // to read the records
    

    public MemorizedIndexedProteinDatabase(String fasta) throws IOException { 
        super(fasta);
        // .pin for peptide index, .min for mass index 
        peptideIndexFile = fasta + ".pin";
        massIndexFile = fasta + ".min";

        timer.startTiming();
        loadMassIndex();
        System.out.println("Time used to load mass index: " + timer.getTimeUsed());        

        MAXINDEX = index.length - 1;
       // b = new byte[DEFAULTNUMRECORDS*RECORDSIZE]; 

        timer.startTiming();
        readPeptideIndex(); 
        System.out.println("Time used to load peptide index: " + timer.getTimeUsed());        
        //raf = new RandomAccessFile(peptideIndexFile, "rws"); 
        //raf.seek(0);
    }
    private void readPeptideIndex() throws IOException {

        File pif = new File(peptideIndexFile); 
        int numpeptides = (int)(pif.length()/RECORDSIZE);
        proteinIndex = new int[numpeptides];    
        peptideStartIndex = new char[numpeptides];    
        peptideEndIndex = new char[numpeptides];    
        System.out.println("Number of Peptide from pin file: " + numpeptides);
        FileInputStream fis = new FileInputStream(pif);
        DataInputStream dis = new DataInputStream(new BufferedInputStream(fis));
        int i = 0;

        System.out.println("dis.available: " + dis.available());
        //while(dis.available() > 0) { // may return negative number if greater than 2G
        while(dis.available() != 0) {
            proteinIndex[i] = dis.readInt();
            peptideStartIndex[i] = dis.readChar();
            peptideEndIndex[i] = dis.readChar();

            i++;
            //System.out.println("dis.available: " + dis.available());
        }

        System.out.println("dis.available: " + dis.available());
        System.out.println("Number of Peptide read from pin file: " + i);
        dis.close();
    } 

    private void diffModSearch(ProcessedPeakList ppl, SearchResult sr, Modifications m) throws IOException {

        // temp solution, will use m.getMassShift()
        //double massShift = m.getDiffNTermMod();
        double massShift = m.getMassShift();

        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();

        timer.startTiming();
        double acc = precMassAccuracy/1000000.0f;
       
        // all static mods except N and C terminal static mods have been 
        // considered while the database are processed, so N and C terminal
        // static mods and all diff mods need special attention here

        double prcMass = ppl.getPrecursorMass();
        
//System.out.println("staticTerminalMods: " + ppl.getSearchParams().getStaticTerminalMods());
        prcMass -= massShift; 
        

        //MassCalculator mc = ppl.getMassCalculator();
        // for multi threading
        //PeptidesReader [] prs = new PeptidesReader[numIsotopes];
        int i = 0; 
        //for(int i = 0; i < numIsotopes; i++) {
        do{
            double mass = 1000*(prcMass - i*MassSpecConstants.MASSDIFFC12C13);
            double diff = mass * acc;
            diff = diff < 500? diff : 500;
            int highLimit = (int)Math.round(mass + diff);
            int lowLimit = (int)Math.round(mass - diff);

            if(numIsotopes == 0) { // traditional sequest like 
                highLimit = (int)Math.round(mass + ppl.getSearchParams().getHighPrecursorTolerance());
                lowLimit = (int)Math.round(mass - ppl.getSearchParams().getLowPrecursorTolerance());
                //System.out.println("highLimit: " + highLimit + "\tlowLimit: " + lowLimit);
            }

            //System.out.print("prcMass: " + prcMass + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);
            long startOffset = index[lowLimit-1];
            if(highLimit > MAXINDEX) { highLimit = MAXINDEX; } 
            long endOffset = index[highLimit];
            int len = (int)(endOffset - startOffset);
//            System.out.println("\tstartOffSet: " + startOffset + "\tendOffset: " + endOffset + "\tlen: " + len);

            if(len > 0) {
                //getPeptides(sr, readPeptides(startOffset, len), m); 
                getPeptides(sr, startOffset/RECORDSIZE, endOffset/RECORDSIZE, m); 
            }
            //prs[i] = new PeptidesReader(sr, startOffset, len); 
            //prs[i].start();
        } while (++i < numIsotopes);
    }
//    public void getPeptideHits(SearchResult sr, MassCalculator mc, double highLimit, double lowLimit) {
    // precMassAccuracy in ppm
    public SearchResult search(ProcessedPeakList ppl) throws IOException {
        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();

        timer.startTiming();
        SearchResult sr = new SearchResult(ppl);
        double acc = precMassAccuracy/1000000.0f;
       
        // all static mods except N and C terminal static mods have been 
        // considered while the database are processed, so N and C terminal
        // static mods and all diff mods need special attention here

        double prcMass = ppl.getPrecursorMass();
        
//System.out.println("staticTerminalMods: " + ppl.getSearchParams().getStaticTerminalMods());
        prcMass -= ppl.getSearchParams().getStaticTerminalMods();
        

        //MassCalculator mc = ppl.getMassCalculator();
        // for multi threading
        //PeptidesReader [] prs = new PeptidesReader[numIsotopes];
        int i = 0; 
        //for(int i = 0; i < numIsotopes; i++) {
        do{
            double mass = 1000*(prcMass - i*MassSpecConstants.MASSDIFFC12C13);
            double diff = mass * acc;
            diff = diff < 500? diff : 500;
            int highLimit = (int)Math.round(mass + diff);
            int lowLimit = (int)Math.round(mass - diff);

            if(numIsotopes == 0) { // traditional sequest like 
                highLimit = (int)Math.round(mass + ppl.getSearchParams().getHighPrecursorTolerance());
                lowLimit = (int)Math.round(mass - ppl.getSearchParams().getLowPrecursorTolerance());
                //System.out.println("highLimit: " + highLimit + "\tlowLimit: " + lowLimit);
            }

           
            if(highLimit > MAXINDEX) { highLimit = MAXINDEX; } 
            if(lowLimit > MAXINDEX) { lowLimit = MAXINDEX; } 
//            System.out.print("highLimit: " + highLimit + "\tlowLimit: " + lowLimit);
            long startOffset = index[lowLimit-1];
            long endOffset = index[highLimit];
            int len = (int)(endOffset - startOffset);
//            System.out.println("\tstartOffSet: " + startOffset + "\tendOffset: " + endOffset + "\tlen: " + len);

            if(len > 0) {
                //getPeptides(sr, readPeptides(startOffset, len)); 
                getPeptides(sr, startOffset/RECORDSIZE, endOffset/RECORDSIZE); 
            }
            //prs[i] = new PeptidesReader(sr, startOffset, len); 
            //prs[i].start();
        } while (++i < numIsotopes);
        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next();
            if(m != null && m.getDiffModsShift() != 0 ) {
//System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                diffModSearch(ppl, sr, m);
            }

        }
//System.out.println("diffmodshift: " + m.getDiffModsShift());
        return sr;
    }

    // for diff mods search
    private void getPeptides(SearchResult sr, long start, long end, Modifications m) throws IOException {
     
        for(int i = (int)start; i <= (int)end; i++) {
            //sr.addPeptideHit(new ModifiedPeptideHit(getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i], m));
            addPeptideHit(sr, getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i], m);
        }       
    }
    private void getPeptides(SearchResult sr, long start, long end) throws IOException {
        for(int i = (int)start; i <= (int)end; i++) {

            //    System.out.println("processing " + proteinIndex[i] + " for " + sr);
            //sr.addPeptideHit(new PeptideHit(getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i]));
            addPeptideHit(sr, getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i]);
        }  
    }
    private void loadMassIndex() throws IOException {
        DataInputStream dis = new DataInputStream(new BufferedInputStream
                                     (new FileInputStream(massIndexFile)));
        int numIndex = dis.available()/8; // 64 is number of byte for a long
        //System.out.println("numIndex: " + numIndex);
        //timer.startTiming();
        index = new long[numIndex];
        for(int i = 0; i < numIndex; i++) {
            index[i] = dis.readLong();
            //if(i > 1587000)
           // System.out.println(i +"\t" + index[i]);
        }
        dis.close();
        //System.out.println("Time used to load mass index: " + timer.getTimeUsed());        
    }
    public static void main(String args[]) throws Exception {
        try { 
            TimeUtils timer = new TimeUtils();
            String fasta = args[0];
            timer.startTiming();
            MemorizedIndexedProteinDatabase se = new MemorizedIndexedProteinDatabase(fasta);
            // peptide index file
            //se.outputPeptides(); 
            //se.sortPeptides(); 
            for(int i = 100000; i < se.index.length; i += 100000) { 
                long size = se.index[i] - se.index[i-1000];
                System.out.println("mass: " + i/1000 + "\tsize: " + size);
            }
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

