/**
 * @file IndexedProteinDatabase.java
 * This is the source file for edu.scripps.pms.mspid.IndexedProteinDatabase
 * @author Tao Xu
 * @date $Date: 2011/11/10 23:42:31 $
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

public class IndexedProteinDatabase extends ProteinDatabase {
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
    private static final String USAGE = "java IndexedProteinDatabase dbfile";
    private String peptideIndexFile;
    private String massIndexFile;
    private RandomAccessFile raf;
    private static final int ACCURACYFACTOR = 1000;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    // for the frequence of peptide length
    private static final int RECORDSIZE = 8;
    //private static final int DEFAULTNUMRECORDS = 1000000; 
    private long [] index;
    private int MAXINDEX;

    long totalNumPeptides = 0;
    // the following variable will be used by the heap sort for swapping
    //private byte [] b = null; // to read the records
    

    public IndexedProteinDatabase(String fasta) throws IOException { 
        super(fasta);
        // .pin for peptide index, .min for mass index 
        peptideIndexFile = fasta + ".pin";
        massIndexFile = fasta + ".min";
        loadMassIndex();
        MAXINDEX = index.length - 1;
       // b = new byte[DEFAULTNUMRECORDS*RECORDSIZE]; 

        raf = new RandomAccessFile(peptideIndexFile, "rws"); 
        raf.seek(0);
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
                getPeptides(sr, readPeptides(startOffset, len), m); 
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
                getPeptides(sr, readPeptides(startOffset, len)); 
            }
            //prs[i] = new PeptidesReader(sr, startOffset, len); 
            //prs[i].start();
        } while (++i < numIsotopes);
        // for multi threading
        /*
        for(int i = 0; i < numIsotopes; i++) {

            try {
               // prs[i].join();
            } catch (Exception e) {
                throw new RuntimeException("Failed to join");
            }

        }
        //System.out.println("Finished join");
        */
        //System.out.println("Time used to get peptide hits: " + timer.getTimeUsed());        
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
    private synchronized byte [] readPeptides(long pos, int len) throws IOException {
        byte [] bytes = new byte[len];
//System.out.println("pos: " + pos + "\tlen: " + len);
        raf.seek(pos);
        raf.read(bytes);
        return bytes;
    }
    // for diff mods search
    private void getPeptides(SearchResult sr, byte [] bytes, Modifications m) throws IOException {

        int numPeptides = bytes.length/RECORDSIZE;
        //System.out.println("numPeptides: " + numPeptides);
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
        for(int i = 0; i < numPeptides; i++) {
            //int mass = dis.readInt();

            //int proteinIndex = dis.readInt();
            //int start = dis.readChar();
            //int end = dis.readChar();  // temporary solution, will remove the minus one 
            //Fasta f = getFasta(proteinIndex);
            //sr.addPeptideHit(new PeptideHit(f, start, end));
            //sr.addPeptideHit(new ModifiedPeptideHit(getFasta(dis.readInt()), dis.readChar(), dis.readChar(), m));
            addPeptideHit(sr, getFasta(dis.readInt()), dis.readChar(), dis.readChar(), m);
        }       
    }
    private void getPeptides(SearchResult sr, byte [] bytes) throws IOException {

        int numPeptides = bytes.length/RECORDSIZE;
        //System.out.println("numPeptides: " + numPeptides);
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
        for(int i = 0; i < numPeptides; i++) {
            //int mass = dis.readInt();

            //int proteinIndex = dis.readInt();
            //int start = dis.readChar();
            //int end = dis.readChar();  // temporary solution, will remove the minus one 
            //Fasta f = getFasta(proteinIndex);
            //sr.addPeptideHit(new PeptideHit(f, start, end));
            //sr.addPeptideHit(new PeptideHit(getFasta(dis.readInt()), dis.readChar(), dis.readChar()));
            addPeptideHit(sr, getFasta(dis.readInt()), dis.readChar(), dis.readChar());
        }       
    }
    private void loadMassIndex() throws IOException {
        DataInputStream dis = new DataInputStream(new BufferedInputStream
                                     (new FileInputStream(massIndexFile)));
        int numIndex = dis.available()/8; // 64 is number of byte for a long
        //System.out.println("numIndex: " + numIndex);
        timer.startTiming();
        index = new long[numIndex];
        for(int i = 0; i < numIndex; i++) {
            index[i] = dis.readLong();
            //if(i > 1587000)
           // System.out.println(i +"\t" + index[i]);
        }
        dis.close();
        System.out.println("Time used to load mass index: " + timer.getTimeUsed());        
    }
    public static void main(String args[]) throws Exception {
        try { 
            TimeUtils timer = new TimeUtils();
            String fasta = args[0];
            timer.startTiming();
            IndexedProteinDatabase se = new IndexedProteinDatabase(fasta);
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

