/**
 * @file ProcessDbNew.java
 * This is the source file for edu.scripps.pms.mspid.ProcessDbNew
 * @author Tao Xu
 * @date $Date: 2011/11/10 23:42:31 $
 */
package edu.scripps.pms.mspid.db;

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.ByteArrayConverter;
import java.util.*;
import java.io.*;
//import java.lang.reflect.Array;

// command line processor
import org.apache.commons.cli.*;


class ProcessDbNew {
    private static final String USAGE= "---Usage: java ProcessDbNew dbfile ---";
    private static final String DBNAMEOPT = "d";
    private static final String TAXONOPT = "t";
    private static final String MOLTYPEOPT = "m";
    private static final String OUTPUTOPT = "o";
    private Protease protease;     
    private int maxMisCleavage = 2;     
    private int enzymespecificity = 0; 
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
    private String peptideIndexFile;
    private String massIndexFile;
    private String fastaFile;
    //private RandomAccessFile raf;
    private static int ACCURACYFACTOR = MassCalculator.MASSACCURACYFACTOR;
    private int AMUSPERFILE = ACCURACYFACTOR/10;
    private int LOWLIMIT; 
    private int HIGHLIMIT;
    private final int LENGTHLIMIT = 6;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    private static int NUMBINS;
    // for the frequence of peptide length
    private int [] massFreq;
    private static final int RECORDSIZE = 12;
    private static final int FINALRECORDSIZE = RECORDSIZE - 4; // ignore mass
//    private byte [] tempNewRecord = new byte[FINALRECORDSIZE];
    private File [] tempFiles;
    long [] massIndex;
    //private static final int NUMRECORDSPERCHUNK = 10000000;
    private static int maxLength = 0;
    private SearchParams sp = null;
    private MassCalculator mc;
    long totalNumPeptides = 0;
    long bytesOutput = 0;
    private PrintStream log;
    //private byte [] temp = new byte[RECORDSIZE];
    // the following variable will be used by the heap sort for swapping
    //private String [] fileNames;
    public ProcessDbNew(String fasta, SearchParams sp, int msicleav) throws IOException { 
        log = new PrintStream(new BufferedOutputStream(new FileOutputStream(fasta + ".processdb.log")));
        this.fastaFile = fasta;
        this.sp = sp;
        maxMisCleavage = msicleav;
        mc = new MassCalculator(sp);
        protease = sp.getProtease();
        enzymespecificity = sp.getEnzymeSpecificity();

        // .pin for peptide index, .min for mass index 
        peptideIndexFile = fastaFile + ".pin";
        massIndexFile = fastaFile + ".min";

    }
    public void printLog() {

        log.println("Database: " + fastaFile);
        log.println("High mass limit: " + HIGHLIMIT);
        log.println("Low mass limit: " + LOWLIMIT);
        log.println("Max num internal miscleavage: " + maxMisCleavage);
        log.println("Peptide index file: " + peptideIndexFile);
        log.println("Mass index file: " + massIndexFile);
        Iterator<Modification> it = sp.getStaticMods();
        if(it.hasNext()) {
            log.println("Static modifications: "); 
        }
        while(it.hasNext()) {
            Modification m = it.next();
            log.println("\t" + m.getResidue() + ": " + m.getMassShift());
        } 
        
    }
    public void printLog(String s) {
        log.println(s);
        
    }
    public void closeLogFile() throws IOException {
        log.close();
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
        int numFiles = (HIGHLIMIT - LOWLIMIT)/(ACCURACYFACTOR*AMUSPERFILE) + 1; 
        String tempdir = "myTemPDir" + System.nanoTime(); 
        File tempDir = new File(tempdir);
        tempDir.deleteOnExit();
        tempDir.mkdir();
        
        tempFiles = new File[numFiles];
        for(int i = 0; i < numFiles; i++) {
            tempFiles[i] = new File(tempDir, ""+ i);
            tempFiles[i].deleteOnExit();
        }
    }

    private void outputMassIndex(long [] massIndex) throws IOException {

        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream
                                     (new FileOutputStream(massIndexFile)));
        for(int i = 0; i < massIndex.length; i++) {
            dos.writeLong(massIndex[i]);
        }
        dos.close();
    }
    private void sortPeptideIndex() throws Exception {

        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream 
                                 (new FileOutputStream(peptideIndexFile)));

        TimeUtils timer = new TimeUtils();
        timer.startTiming();
        int lastMass = 0;
        for(File f : tempFiles) {
            //System.out.println("Now sorting " + f);
            lastMass = sortFile(f, dos, lastMass);
        }
        outputMassIndex(massIndex);
        timer.stopTiming();
        System.out.println("Time used for sorting: " + timer.getTimeUsed());
        //System.out.println("After sorting");
        System.gc();
    }
    private void outputPeptides() throws IOException {
        TimeUtils timer = new TimeUtils();
        timer.startTiming();
        DataOutputStream [] dos = new DataOutputStream[tempFiles.length];
        for(int i = 0; i < tempFiles.length; i++) {
            dos[i] = new DataOutputStream(new BufferedOutputStream
                                     (new FileOutputStream(tempFiles[i])));
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
        printLog("Number of proteins in the database: " + i + "\tLongest protein: " + maxLength);
        printLog("Number of peptides: " + totalNumPeptides + "\tavgNumPeptidesPerProtein: " + totalNumPeptides/(double)i);
        printLog("NumZeros: " + numZeros + "\tNumNonZeros: " + numNonZeros);
        fis.close();
        for(i = 0; i < dos.length; i++) {
            dos[i].close();
        }
    }

    
    // index is the index for the fasta in the database
    private int outputPeptides(int index, Fasta f, DataOutputStream [] dos) throws IOException {
        
         // to check the longest protein in the database
        int lastIndex = f.getLength();
        if(lastIndex > maxLength) { 
            maxLength = lastIndex;
        }
        byte [] seq = f.getSequenceAsBytes();
        int numPeptides = 0; 
        for(char i = 0; i < lastIndex; i++) {
            char j = i; // use char as unsigned short
             // mass of H2O need to be added
             // because the mass of the residues does not count H2O
            float mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSPROTON;
            int length = 0;
            while(mass*ACCURACYFACTOR <= HIGHLIMIT && j < lastIndex) {
                
                mass += mc.getPrecursorMass(seq[j]);
                length++;
                int intMass = Math.round(mass * ACCURACYFACTOR);
                if(length >= LENGTHLIMIT && intMass >= LOWLIMIT && intMass <= HIGHLIMIT) {
                    //if(protease.checkEnzymeSpecificityStrict(seq, i, j) >= enzymespecificity) {
                    if(protease.checkEnzymeSpecificityStrict(seq, i, j) >= enzymespecificity && 
                       //protease.getNumInternalMissCleavage(seq, i, j) <= 2 ) {
                       protease.getNumInternalMissCleavage(seq, i, j) <= maxMisCleavage ) {
                        int dosindex = (intMass-LOWLIMIT)/(AMUSPERFILE*ACCURACYFACTOR);
                        dos[dosindex].writeInt(intMass); // mass
                        dos[dosindex].writeInt(index);  // index of the fasta
                        dos[dosindex].writeChar(i);  // start of the peptide
                        dos[dosindex].writeChar(j); // end of the peptide
                        numPeptides++;
                        massFreq[intMass]++; 
                    }
                }
                j++;
            }
        }
        return numPeptides;
    }

    // not used
    private static CommandLine parseCommandLine(String args[]) {

        // Command line processing
        Options opts = new Options();
        Option dbNameOpt = new Option 
            (DBNAMEOPT, "database", true, 
             "Name of the fasta file. e.g., refseqr4.");
        Option taxonOpt = new Option
            (TAXONOPT, "taxonid", true,
             "An int for taxonomy id. e.g., 9606 for human");
        Option moltypeOpt = new Option
            (MOLTYPEOPT, "moltype", true, 
             "An int for molecular type. e.g., 3 for protein.");
        Option outputOpt = new Option
            (OUTPUTOPT, "output", true, 
             "The name of the output file.");
        dbNameOpt.setRequired(true);
        taxonOpt.setRequired(true);
        moltypeOpt.setRequired(true);
        outputOpt.setRequired(true);
        opts.addOption(dbNameOpt);
        opts.addOption(taxonOpt);
        opts.addOption(moltypeOpt);
        opts.addOption(outputOpt);

        BasicParser cliParser = new BasicParser();
        CommandLine cli = null;

        try {
            cli = cliParser.parse(opts, args);
        } catch (Exception e) {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp(USAGE, opts);

            System.out.println
                ("\nNB: Please use \"java -Xmx500M ... \".");
            System.exit(1);
        }

        return cli;

    }

    public static void main(String args[]) throws Exception {
        try { 
            TimeUtils timer = new TimeUtils();
            SearchParams sp = new SearchParams("search.xml");
            String fastaFile = args[0];
            int maxmiscleavage = 1;
            if(args.length > 1) {
                maxmiscleavage = Integer.parseInt(args[1]);
            }
            timer.startTiming();
            ProcessDbNew se = new ProcessDbNew(fastaFile, sp, maxmiscleavage);
            se.setParameters(args);
            se.printLog();
            se.outputPeptides();
            se.sortPeptideIndex();
            // prepare for garbage collection to reduce memory usage
            System.out.println("Finished sorting");
            //se.writeMassIndex();

            timer.stopTiming();
            long timeUsed = timer.getTimeUsedMillis(); 
            System.out.println("Time used for process database: " + timeUsed);
            se.printLog("Time used for process database: " + timeUsed);
            se.closeLogFile();
            System.out.println();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
        /**
         * @return return true if the temp array is aready sorted before calling this function
         */
        private int sortFile(File f, DataOutputStream dos, int lastMass) throws IOException {
            TimeUtils timer = new TimeUtils(); timer.startTiming();
            int numPeptides = (int)f.length()/RECORDSIZE;
            if(numPeptides < 1) { return lastMass; }
            DataInputStream dis = new DataInputStream(new BufferedInputStream
                                           (new FileInputStream(f)));
            ArrayList<Peptide> peptides = new ArrayList<Peptide>(numPeptides);
            for(int i = 0; i < numPeptides; i++) {
                peptides.add(new Peptide(dis.readInt(), dis.readInt(), dis.readChar(), dis.readChar()));
            } 
            dis.close();
            timer.stopTiming();
            System.out.println("Time used to read peptides in file " + f + ": " + timer.getTimeUsed());
            timer.startTiming();

            Collections.sort(peptides);
            timer.stopTiming();
            System.out.println("Time used to sort peptides: " + timer.getTimeUsed());
            timer.startTiming();
            int currentMass = 0;            
            for(Peptide p : peptides) {
                //dos.writeInt(p.mass);
                dos.writeInt(p.proteinIndex);
                dos.writeChar(p.start);
                dos.writeChar(p.stop);
                bytesOutput += FINALRECORDSIZE;
                if(p.mass > lastMass) {
                    //long index = dos.size() - FINALRECORDSIZE;
                    long index = bytesOutput - FINALRECORDSIZE;
                    for(int i = lastMass; i < p.mass; i++) {
                        massIndex[i] = index; 
                    }
                    lastMass = p.mass;
                }
                currentMass = p.mass;
            }
            for(int i = lastMass; i <= currentMass; i++) {
                massIndex[i] = bytesOutput;
            } 
            peptides = null; System.gc();
            timer.stopTiming();
            System.out.println("Time used to output peptide index: " + timer.getTimeUsed());
            return currentMass;
        }
}


