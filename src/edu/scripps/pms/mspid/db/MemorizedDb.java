/**
 * @file MemorizedDb.java
 * This is the source file for edu.scripps.pms.mspid.MemorizedDb
 * @author Tao Xu
 * @date $Date: 2011/11/10 23:42:31 $
 */
package edu.scripps.pms.mspid.db;

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.ByteArrayConverter;
//import edu.scripps.pms.util.enzyme.Protease;
import java.util.*;
import java.io.*;
import java.lang.reflect.Array;

public class MemorizedDb extends ProteinDatabase {
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
    private static final String USAGE = "java MemorizedDb dbfile";
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
    private SearchParams params = null;
    private Protease protease;     
    private int maxMisCleavage = 2; // need to read from parameter file 
    private MassCalculator mc;
    private double maxMass;
    private double minMass;
  

    public MemorizedDb(String fasta, SearchParams sp, ArrayList<Double> prcmasses) throws IOException { 
        super(fasta);

         maxMass = Collections.max(prcmasses); // need to consider diffmod
         minMass = Collections.min(prcmasses);

        params = sp;
        
        //maxMisCleavage = msicleav;
        mc = new MassCalculator(sp);
        protease = sp.getProtease();
        indexDatabase(prcmasses);

    }

    private void processMass(boolean [] bmasses, double mass) {
        double acc = params.getPrecursorTolerance()/1000000.0f;        
        double diff = mass * acc;
      
System.out.println(mass);
        for(int i = 0; i < params.getNumIsotopicPeaks(); i++) {
            double tempmass = mass - i*MassSpecConstants.MASSDIFFC12C13;
            int lowintmass = (int)((tempmass - diff)*1000);
            int highintmass = (int)((tempmass + diff)*1000);
//System.out.println(lowintmass + "\t" + highintmass + "\t" + diff + "\t" + acc + "\t" + mass);
             for(int intmass = lowintmass; intmass <= highintmass; intmass++) {
  
                 bmasses[intmass] = true;
             }
        }
    }

    private boolean [] getInterestedMasses(ArrayList<Double> prcmasses) {
        boolean [] interested = new boolean[100000000]; // handles 100,000 D
        Iterator<Double> it = prcmasses.iterator();        
        while(it.hasNext()){
            double prcmass = it.next();
            processMass(interested, prcmass);
            
            //Iterator<DiffMod> diffmodit = params.getDiffMods();
            Iterator<Modifications> diffmodit = params.getAllModifications();
            while(diffmodit.hasNext()) {
                Modifications dm = diffmodit.next();
                double massshift = dm.getDiffModsShift();
                if(massshift != 0) {
                    double mass = prcmass - dm.getMassShift();
                
                    processMass(interested, mass);
                }
            }

        }
        return interested;
    }
    public void indexDatabase(ArrayList<Double> prcmasses) {
System.out.println("Total Number of precursors: " + prcmasses.size());
        boolean [] bmasses = getInterestedMasses(prcmasses);
        int numProteins = 0;
        int totalNumPeptides = 0;
        ArrayList<edu.scripps.pms.mspid.Peptide> [] peptides = new ArrayList[10000000];
        int enzymespecificity = params.getEnzymeSpecificity();
        int minPeptideLength = params.getMinimumPeptideLength();
        Iterator<Fasta> proteins = sequences.iterator();
        int [] massFreq = new int[10000000];
        while(proteins.hasNext()) {
            Fasta f = proteins.next();
            numProteins++;


            int lastIndex = f.getLength();
            byte [] seq = f.getSequenceAsBytes();
            int numPeptides = 0; 

            for(int i = 0; i < lastIndex; i++) {
                int j = i; // use char as unsigned short
                // mass of H2O need to be added
                // because the mass of the residues does not count H2O
                // need to consider n-term and c-term static mods
                double mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSPROTON + params.getStaticTerminalMods();
                int length = 0;
                while(mass <= maxMass && j < lastIndex) {
                    
                    mass += mc.getPrecursorMass(seq[j]);
                    length++;
                    int intMass = (int)(mass * ACCURACYFACTOR);
                    if(bmasses[intMass]) { 
                        if(length >= minPeptideLength && mass >= minMass && mass <= maxMass) {
                            //if(protease.checkEnzymeSpecificityStrict(seq, i, j) >= enzymespecificity) {
                            if(protease.checkEnzymeSpecificityStrict(seq, i, j) >= enzymespecificity && 
                                                protease.getNumInternalMissCleavage(seq, i, j) <= 2 ) {
                                numPeptides++;
                                massFreq[intMass]++; 
                                if(peptides[intMass] == null) {
                                    peptides[intMass] = new ArrayList<edu.scripps.pms.mspid.Peptide>();
                                }
                                edu.scripps.pms.mspid.Peptide p = new edu.scripps.pms.mspid.Peptide(f, i, j);
                                peptides[intMass].add(p);
                                numPeptides++;
                            }
                        }
                    }
                    j++;
                }
            }

            totalNumPeptides += numPeptides;

            System.out.println("Finished processing " + f.getAccession() + "\tNumProtiens " + numProteins + "\tNumPeptides " + numPeptides + "\tTotalNumPeptides " + totalNumPeptides); 

 
        }
 
        System.out.println("Finished processing database indexing. NumProtiens processed: " + numProteins + "\tTotalNumPeptides: " + totalNumPeptides); 
        int maxfreq = 0;
        for(int i = 0; i < massFreq.length; i++) {
            if(massFreq[i] > maxfreq) {
                maxfreq = massFreq[i];
            }
        }
        System.out.println("MaxFreq: " + maxfreq);

    }

    private void diffModSearch(ProcessedPeakList ppl, SearchResult sr, Modifications m) throws IOException {

        // temp solution, will use m.getMassShift()
        //double massShift = m.getDiffNTermMod();
        double massShift = m.getMassShift();

        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();

        timer.startTiming();
        double acc = precMassAccuracy/1000000.0f;
       
        // all static mods except N and C terminal static mods have been 
        // considered while the database are processed, so N and C terminal
        // static mods and all diff mods need special attention here

        double prcMass = ppl.getPrecursorMass();
        
        prcMass -= massShift; 
        

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
        
        prcMass -= ppl.getSearchParams().getStaticTerminalMods();
        

        int i = 0; 
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
            //sr.addPeptideHit(new ModifiedPeptideHit(getFasta(dis.readInt()), dis.readChar(), dis.readChar(), m));
            addPeptideHit(sr, getFasta(dis.readInt()), dis.readChar(), dis.readChar(), m);
        }       
    }
    private void getPeptides(SearchResult sr, byte [] bytes) throws IOException {

        int numPeptides = bytes.length/RECORDSIZE;
        //System.out.println("numPeptides: " + numPeptides);
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
        for(int i = 0; i < numPeptides; i++) {
            //sr.addPeptideHit(new PeptideHit(getFasta(dis.readInt()), dis.readChar(), dis.readChar()));
            addPeptideHit(sr, getFasta(dis.readInt()), dis.readChar(), dis.readChar());
        }       
    }
    public static void main(String args[]) throws Exception {
        try { 
            TimeUtils timer = new TimeUtils();
            timer.startTiming();
            String searchxml = "search.xml";
            SearchParams sp = new SearchParams(searchxml);
            String fasta = sp.getDbName();
            String msFileName = ".";
            ArrayList<String> ms2files = ProlucidSearchEngine.getMs2Files(msFileName);

            ArrayList<Double> allprecursormasses = ProlucidSearchEngine.getAllPrecursorMasses(ms2files, sp);
 
            MemorizedDb se = new MemorizedDb(fasta, sp, allprecursormasses);
            long totalmemory = Runtime.getRuntime().totalMemory();
            long usedmemory  = totalmemory - Runtime.getRuntime().freeMemory();
            
            System.out.println("Total memory\t" + totalmemory + "\n");
            System.out.println("Memory used after loading the database\t" + usedmemory + "\n");
            // peptide index file
            //se.outputPeptides(); 
            //se.sortPeptides(); 
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

