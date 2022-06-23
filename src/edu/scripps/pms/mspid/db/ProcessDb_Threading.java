/**
 * @file ProcessDb_Thread.java
 * This is the source file for edu.scripps.pms.mspid.ProcessDb_Thread
 * @author Tao Xu
 * @date $Date: 2010/11/04 00:06:06 $
 */
package edu.scripps.pms.mspid.db;

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.ByteArrayConverter;
import java.util.*;
import java.io.*;
//import java.lang.reflect.Array;

// command line processor
import org.apache.commons.cli.*;


class ProcessDb_Thread {
    private static final String USAGE= "---Usage: java ProcessDb_Thread dbfile ---";
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
    private String peptideIndexFile;
    private String massIndexFile;
    private String fastaFile;
    //private RandomAccessFile raf;
    private static int ACCURACYFACTOR = 1000;
    private int AMUSPERFILE = 100;
    private int LOWLIMIT; 
    private int HIGHLIMIT;
    private final int LENGTHLIMIT = 7;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    private static int NUMBINS;
    // for the frequence of peptide length
    private int [] massFreq;
    private static final int RECORDSIZE = 12;
    private static final int FINALRECORDSIZE = RECORDSIZE - 4; // ignore mass
//    private byte [] tempNewRecord = new byte[FINALRECORDSIZE];
    private File [] tempFiles;
    //private static final int NUMRECORDSPERCHUNK = 10000000;
    private static int maxLength = 0;
    private SearchParams sp = null;
    private MassCalculator mc;
    long totalNumPeptides = 0;
    // the following variable will be used by the heap sort for swapping
    //private String [] fileNames;
    public ProcessDb_Thread(String fasta, SearchParams sp) throws IOException { 
        this.fastaFile = fasta;
        this.sp = sp;
        mc = new MassCalculator(sp);
        // .pin for peptide index, .min for mass index 
        peptideIndexFile = fastaFile.substring(0, fastaFile.lastIndexOf(".")) + ".pin";
        massIndexFile = fastaFile.substring(0, fastaFile.lastIndexOf(".")) + ".min";

    }
    private void setParameters(String [] args) {

        LOWLIMIT = (int)sp.getMinPrecursorMass()*ACCURACYFACTOR;
        HIGHLIMIT = (int)sp.getMaxPrecursorMass()*ACCURACYFACTOR;
        // keep 1pp accuracy for peptides with mass 1000 dalton
        NUMBINS = HIGHLIMIT + 1; 
        // for the frequence of peptide length
        massFreq = new int[NUMBINS];
        // one file for every AMUSPERFILE amu
        int numFiles = (HIGHLIMIT - LOWLIMIT)/(ACCURACYFACTOR*AMUSPERFILE) + 1; 
         
        tempFiles = new File[numFiles];
        for(int i = 0; i < numFiles; i++) {
            tempFiles[i] = new File("temp" + i);
            tempFiles[i].deleteOnExit();
        }
    }
    private void writeMassIndex() throws IOException {

        long [] massIndex = new long[NUMBINS];
        int lastMass = 0;
        int currentMass = 0;
        long bytesRead = 0;      
        TimeUtils timer = new TimeUtils(); 
        timer.startTiming();
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream 
                                 (new FileOutputStream(peptideIndexFile)));
        for(File f : tempFiles) {
            // to avoid creating the temp array  too many times
            byte [] temp = new byte[FINALRECORDSIZE]; 
            DataInputStream dis = new DataInputStream(new BufferedInputStream 
                                 (new FileInputStream(f)));
            int numBytes = dis.available();
            if(numBytes == 0) continue;
//System.out.println("Processing file temp" + f + "\tnumBytes: " + numBytes);
            
            int numRecordsProcessed = 0;
            
            int numPeptides = numBytes/RECORDSIZE;
            while(numRecordsProcessed < numPeptides) {
                while(currentMass == lastMass && numRecordsProcessed < numPeptides) {

//System.out.println("numPeptides: " + numPeptides + "\tnumRecrodsProcessed: " + numRecordsProcessed);
                    currentMass = processRecord(dis, dos, temp); 
                    numRecordsProcessed++;
                    bytesRead += FINALRECORDSIZE;
                }
                //System.out.println("CurentMass: " + currentMass + "\tlastMass: " + lastMass);
                if(currentMass < lastMass) {
    
                    System.out.println("found!, currentMass is: " + currentMass + "\tlastMass: " + lastMass); 
                    throw new RuntimeException("The peptide index file is not sorted by mass");
                }               
                for(int k = lastMass; k < currentMass; k++) {
                    massIndex[k] = bytesRead - FINALRECORDSIZE;

                }
                    //bytesRead += RECORDSIZE;
                lastMass = currentMass;
            }
            dis.close();

        }
        for(int i = lastMass; i <= currentMass; i++) {
            massIndex[i] = bytesRead;
        }
        dos.close();
        outputMassIndex(massIndex);
        timer.stopTiming();
        long timeUsed = timer.getTimeUsedMillis(); // end of the last mass
        System.out.println("Time used to writeMassIndex: " + timeUsed);
    }




    private int processRecord(DataInputStream dis, DataOutputStream dos, byte [] temp) throws IOException {

        int mass = dis.readInt();
        dis.read(temp); 
        dos.write(temp);
        //dis.skip(8);
        //dis.readInt(); dis.readChar(); dis.readChar();
        return mass;
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


        //TimeUtils timer = new TimeUtils();
        Thread mainThread = Thread.currentThread();
        int numFilesProcessed = 0;
        while(numFilesProcessed < tempFiles.length) {
            
            if(Thread.activeCount() < 5) {
                new SortingThread(tempFiles[numFilesProcessed++]).start();
            } else {
                mainThread.yield();
            }
        }

        while(Thread.activeCount() > 1) {
            mainThread.yield();
        }
        System.out.println("After sorting");
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
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            totalNumPeptides += outputPeptides(i++, f, dos);
        //    sequences.add(f);
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
                    int dosindex = (intMass-LOWLIMIT)/(AMUSPERFILE*ACCURACYFACTOR);
                    dos[dosindex].writeInt(intMass); // mass
                    dos[dosindex].writeInt(index);  // index of the fasta
                    dos[dosindex].writeChar(i);  // start of the peptide
                    dos[dosindex].writeChar(j); // end of the peptide
                    numPeptides++;
                    massFreq[intMass]++; 
                }
                j++;
            }
        }
        return numPeptides;
    }

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
            timer.startTiming();
            ProcessDb_Thread se = new ProcessDb_Thread(fastaFile, sp);
            se.setParameters(args);

            se.outputPeptides();
            se.sortPeptideIndex();
            // prepare for garbage collection to reduce memory usage
            System.out.println("Finished sorting");
            se.writeMassIndex();

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
    private class SortingThread extends Thread {
        private File file;
        byte [] temp = new byte[RECORDSIZE];
        public SortingThread(File f) {
            file = f;
        }
        public void run() {
            try {
                System.out.println("sorting file " + file.getName());

                TimeUtils timer = new TimeUtils();
                timer.startTiming();
                sortFile(file);
                timer.stopTiming();
                long timeUsed = timer.getTimeUsed();
                System.out.println("TimeUsed: " + timeUsed + "\tfor sorting file: "+ file.getName());
                System.gc();
                System.out.println("finished sorting file " + file.getName());
                } catch (IOException e) {
                    e.printStackTrace();
                    System.exit(1);
                }
        }
        /**
         * @return return true if the temp array is aready sorted before calling this function
         */
        private boolean sortFile(File f) throws IOException {
            int numPeptides = (int)f.length()/RECORDSIZE;
System.out.println("File length before write: " + f.length());
            DataInputStream dis = new DataInputStream(new BufferedInputStream
                                           (new FileInputStream(f)));
            
            ArrayList<Peptide> peptides = new ArrayList<Peptide>(numPeptides);
            for(int i = 0; i < numPeptides; i++) {
                int mass = dis.readInt();
//System.out.println("mass: " + mass);
                byte [] location = new byte[FINALRECORDSIZE];
                dis.read(location);
                peptides.add(new Peptide(mass, location));
            } 
            dis.close();
System.out.println("File length after read: " + f.length());
            if(peptides.size() > 0) {
System.out.println("numpeptides added: " + numPeptides);
System.out.println("numpeptides write: " + peptides.size());
                Collections.sort(peptides);
//System.out.println("First mass: " + peptides.get(0).mass + "\tlast mass: " + peptides.get(peptides.size()-1).mass);
                DataOutputStream dos = new DataOutputStream(new BufferedOutputStream
                                          (new FileOutputStream(f)));
                for(Peptide p : peptides) {
                    dos.writeInt(p.mass);
                    dos.writeInt(p.proteinIndex);
                    dos.writeChar(p.start);
                    dos.writeChar(p.stop);
                }
                dos.close();
            }
System.out.println("File length after write: " + f.length());
            return false;

/*
            RandomAccessFile raf = new RandomAccessFile(f, "rws");
            int length = (int) raf.length();
            byte [] bytes = new byte[length];
            raf.read(bytes);
            if(checkSorted(bytes)) {
                return true;
            }
            sortRecords(bytes);
            if(!checkSorted(bytes)) {
                System.out.println("!!!!!!!! Not correctly sorted !!!!!!!!!");
                // throw new RuntimeException("!!!!!!!! Not correctly sorted !!!!!!!!!");
            } else {
                System.out.println("Correctly sorted");
            }
            raf.seek(0);
            raf.write(bytes);
            raf.close();

            raf.close();
            raf = null;
            bytes = null;
            System.gc();
            return false;
*/
        }
        /**
         * Heap sort algorithm is used
         * @param b - the byte array to be sorted
         * @param recordSize - the size of each record
         * @param startByte - the index of the first byte of the field to sort upon 
         * @param numBytes - number of bytes of the field to sort upon. e.g., 4 for int
         */
        private void sortRecords(byte [] bytes) {
            int numRecords = bytes.length/RECORDSIZE;
            for(int i = (numRecords/2)-1; i >= 0; i--) {
                siftDown(i, numRecords-1, bytes);
            }
    
            for(int i = numRecords-1; i >= 1; i--) {
                swap(0, i, bytes);
                siftDown(0, i-1, bytes);
            }
        }
        // replace the contents of the record stored in a by the content in b
        // temp is an array of 4 bytes 
        private void copyIn(int index, byte [] bytes) {
            index *= RECORDSIZE;
            for(byte i = 0; i < RECORDSIZE; i++) {
                bytes[index++] = temp[i];
            }
        }
        private void copyOut(int index, byte [] bytes) {
            index *= RECORDSIZE;
            for(byte i = 0; i < RECORDSIZE; i++) {
                temp[i] = bytes[index++];
            }
        }
        private void swap(int indexa, int indexb, byte [] bytes) {
            copyOut(indexa, bytes); // copy data out into temp
            int ia = indexa*RECORDSIZE;
            int ib = indexb*RECORDSIZE;
            for(byte i = 0; i < RECORDSIZE; i++) {
                bytes[ia++] = bytes[ib++];
            }
            copyIn(indexb, bytes);  // copy data in from temp
        }

        // return true if it is already sorted
        private void siftDown(int root, int bottom, byte [] bytes) {
            boolean done = false;
            int maxChild = 0;
            while((root*2 <= bottom) && !done) {
                if(root*2 == bottom) {
                    maxChild = root * 2;
                } else if(compare(root*2, root*2+1, bytes) > 0) {
                    maxChild = root * 2;
                } else {
                    maxChild = root * 2 + 1;
                }
                if (compare(root, maxChild, bytes) < 0) {
                    // swap
                   swap(root, maxChild, bytes);
               
                    // update root 
                    root = maxChild;
                } else {
                    done = true;
                }

            }
        }
        private boolean checkSorted(byte [] bytes) {
            for(int i = 1; i < bytes.length/RECORDSIZE; i++) {
                if(compare(i, i-1, bytes) < 0) {
                    return false;
                }
            }
            return true;
        }
        private int compare(int indexA, int indexB, byte [] bytes) {
            int ia = indexA * RECORDSIZE;
            int ib = indexB * RECORDSIZE;

            if(bytes[ia] != bytes[ib]) return (bytes[ia]&0xff) - (bytes[ib]&0xff);
            if(bytes[++ia] != bytes[++ib]) return (bytes[ia]&0xff) - (bytes[ib]&0xff);
            if(bytes[++ia] != bytes[++ib]) return (bytes[ia]&0xff) - (bytes[ib]&0xff);
            if(bytes[++ia] != bytes[++ib]) return (bytes[ia]&0xff) - (bytes[ib]&0xff);

            return 0;
        }
    }
}


