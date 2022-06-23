/**
 * @file ProcessDb.java
 * This is the source file for edu.scripps.pms.mspid.ProcessDb
 * @author Tao Xu
 * @date $Date: 2009/02/02 05:11:43 $
 */
package edu.scripps.pms.mspid.db;

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
//import edu.scripps.pms.util.ByteArrayConverter;
import java.util.Iterator;
import java.io.*;
//import java.lang.reflect.Array;

// command line processor
import org.apache.commons.cli.*;


class ProcessDb {
    private static final String USAGE= "---Usage: java ProcessDb dbfile ---";
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
    private TimeUtils timer = new TimeUtils();
    private static long startTime;
    private static long endTime;
    private String peptideIndexFile;
    private String massIndexFile;
    private String fastaFile;
    //private RandomAccessFile raf;
    private static int ACCURACYFACTOR = 1000;
    private static int LOWLIMIT; 
    private static int HIGHLIMIT;
    private static final int LENGTHLIMIT = 7;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    private static int NUMBINS;
    // for the frequence of peptide length
    private int [] massFreq;
    private static final int RECORDSIZE = 12;
    //private static final int NUMRECORDSPERCHUNK = 10000000;
    private static int NUMRECORDSPERCHUNK = 8388608*2*2;
    private static int NUMBYTESPERCHUNK = RECORDSIZE*NUMRECORDSPERCHUNK;
    private static int NUMBYTESHALFCHUNK = NUMBYTESPERCHUNK/2;
    private static int maxLength = 0;
    //private SearchParams sp = null;
    private MassCalculator mc;
    long totalNumPeptides = 0;
    // the following variable will be used by the heap sort for swapping
    private byte [] temp = new byte[RECORDSIZE];
    private byte [] b = null;
    private byte [] sortedChunk = null;

    public ProcessDb(String fasta, SearchParams sp) throws IOException { 
        this.fastaFile = fasta;
        //this.sp = sp;
        mc = new MassCalculator(sp);
        // .pin for peptide index, .min for mass index 
        peptideIndexFile = fastaFile.substring(0, fastaFile.lastIndexOf(".")) + ".pin";
        massIndexFile = fastaFile.substring(0, fastaFile.lastIndexOf(".")) + ".min";


    }
    private void setParameters(String [] args) {

        LOWLIMIT = 500*ACCURACYFACTOR;
        HIGHLIMIT = 4500*ACCURACYFACTOR;
        // keep 1pp accuracy for peptides with mass 1000 dalton
        NUMBINS = HIGHLIMIT + 1; 
        // for the frequence of peptide length
        massFreq = new int[NUMBINS];

    }
    private long getNumPeptides() throws IOException {

        RandomAccessFile raf = new RandomAccessFile(peptideIndexFile, "rws"); 
        long peptideIndexFileLength = raf.length();
        raf.close();
        if(peptideIndexFileLength%RECORDSIZE != 0) {
            throw new RuntimeException("numBytes must be multiples of recordSize");
        }
        return peptideIndexFileLength/RECORDSIZE;
     
    }
    private void writeMassIndex() throws IOException {
        long numPeptides = getNumPeptides();

        long [] massIndex = new long[NUMBINS];
        int lastMass = 0;
        int currentMass = 0;
        long bytesRead = 0;      
        
        timer.startTiming();
        DataInputStream dis = new DataInputStream(new BufferedInputStream 
                                 (new FileInputStream(peptideIndexFile), NUMBYTESPERCHUNK));
        int numRecordsProcessed = 0;
        while(numRecordsProcessed < numPeptides) {
            while(currentMass == lastMass && numRecordsProcessed < numPeptides) {
                currentMass = readRecord(dis); 
                numRecordsProcessed++;
                bytesRead += RECORDSIZE;
            }
            long bytesProcessed = bytesRead - RECORDSIZE;
            //System.out.println("CurentMass: " + currentMass + "\tlastMass: " + lastMass);
            if(currentMass < lastMass) {
    
                System.out.println("found!, currentMass is: " + currentMass + "\tlastMass: " + lastMass +  "\tbytesProcessed: " + bytesProcessed);
                throw new RuntimeException("The peptide index file is not sorted by mass");
            }               
            for(int k = lastMass; k < currentMass; k++) {
                massIndex[k] = bytesProcessed;
            }
                    //bytesRead += RECORDSIZE;
            lastMass = currentMass;
        }
        for(int i = lastMass; i <= currentMass; i++) {
            massIndex[i] = bytesRead;
        }
        dis.close();
        outputMassIndex(massIndex);
        timer.stopTiming();
        long timeUsed = timer.getTimeUsedMillis(); // end of the last mass
        System.out.println("Time used to writeMassIndex: " + timeUsed);
    }




    private int readRecord(DataInputStream dis) throws IOException {
        int mass = dis.readInt();
        dis.skip(8);
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
    private void setChunkSize(long peptideIndexFileLength) {
        if(peptideIndexFileLength%RECORDSIZE != 0) {
            throw new RuntimeException("numBytes must be multiples of recordSize");
        }
        if(peptideIndexFileLength <= NUMBYTESPERCHUNK) {
            NUMBYTESPERCHUNK = (int)peptideIndexFileLength;
            NUMRECORDSPERCHUNK = NUMBYTESPERCHUNK/RECORDSIZE;
            NUMBYTESHALFCHUNK = NUMBYTESPERCHUNK/2;
        }
        b = new byte[NUMBYTESPERCHUNK];
        sortedChunk = new byte[NUMBYTESPERCHUNK]; 
    }
    private void sortPeptideIndex() throws IOException {
        RandomAccessFile raf = new RandomAccessFile(peptideIndexFile, "rws"); 
        long peptideIndexFileLength = raf.length();
        setChunkSize(peptideIndexFileLength);
        long numPeptides = peptideIndexFileLength/RECORDSIZE;
        long numChunks = numPeptides/NUMRECORDSPERCHUNK;
        long remainder = peptideIndexFileLength%NUMBYTESPERCHUNK;
        System.out.println("NUMRECORDSPERCHUNK: " + NUMRECORDSPERCHUNK + "\tnumChunks: " + numChunks + "\tRemainder: " + remainder);
        //long timeUsed = 0;
        long lastChunkIndex = peptideIndexFileLength - NUMBYTESPERCHUNK;        

        int numRounds = 0; 
        boolean isSorted = false;       
        TimeUtils timer = new TimeUtils();
         
        if(remainder != 0) {
            isSorted &= sortChunk(raf, lastChunkIndex, b);
        }
        for(long i = 0; i < numChunks; i++) {
            long pos = i * NUMBYTESPERCHUNK;
            timer.startTiming();
            isSorted &= sortChunk(raf, pos, b);
            timer.stopTiming();
            long timeUsed = timer.getTimeUsed();
            System.out.println("pos: " + pos + "\ttimeUsed: " + timeUsed + "\tfor sorting chunk number: "+ (i+1) + "\r");
        }
        if(isSorted) {
            System.out.println("Is sorted already");
        } else {
            System.out.println("Not sorted yet");
        }
        System.out.println("After sorting, begin merging");
        while(!isSorted) {
            isSorted = true;
            System.out.println("numRounds: " + ++numRounds + "\r");
            long numShiftChunks = (remainder < NUMBYTESPERCHUNK/2) ? numChunks-1 : numChunks;
            System.out.println("numChunks: " + numChunks + "\tnumShiftChunks: " + numShiftChunks);
            for(long i = 0; i < numShiftChunks; i++) {
                long pos = i * NUMBYTESPERCHUNK + NUMBYTESHALFCHUNK;
                isSorted &= mergeSortedChunks(raf, pos);
            }
            if(remainder != 0) {
                isSorted &= mergeSortedChunks(raf, lastChunkIndex);
            }
            for(long i = 0; i < numChunks; i++) {
                long pos = i * NUMBYTESPERCHUNK;
                isSorted &= mergeSortedChunks(raf, pos);
            }
        }
        b = null;
        sortedChunk = null;
        System.gc();
        raf.close();
    }
    // find the index of the first element that is smaller than the one before it
    // return 0 if the array is aready in proper order
    private int getSwitchIndex(byte [] bytes) {
        int index = 0;
        for(int i = 1; i < NUMRECORDSPERCHUNK; i++) {
            if(compare(i, i-1, bytes) < 0) {
                return i;
            }
        } 
        return index;
    }
        
    /**
     * @return return true if the temp array is aready sorted before calling this function
     * @note this function requre the first and second half are sorted repectively
     */
    private boolean mergeSortedChunks(RandomAccessFile raf, long pos) throws IOException {
        TimeUtils timer = new TimeUtils();
        timer.startTiming();
        raf.seek(pos);
        raf.read(b);
        if(checkSorted(b)) {
    
            System.out.println("Already sorted before merge, time used to check sorted: " + timer.getTimeUsed());
            return true;
        }
        int ff = 0;  // index of first element in first half
        //int fs = NUMRECORDSPERCHUNK/2; // index of firstElement in second half
        int fs = getSwitchIndex(b); // index of firstElement in second half
        if(fs == 0) {
            System.out.println("Already sorted before merge");
            return true; // is sorted already
        } 
        int lf = fs - 1; // last element in first half
        int ls = NUMRECORDSPERCHUNK - 1; // last element in the last half
        int numRecordMerged = 0;
        while (numRecordMerged < NUMRECORDSPERCHUNK) {
            if(ff <= lf && fs <= ls) {
                if(compare(ff, fs, b) < 0) { // ff is not greater than fs
                    copyRecord(ff++, numRecordMerged++, b, sortedChunk);
                } else {
                    copyRecord(fs++, numRecordMerged++, b, sortedChunk);
                }
            } else if(ff <= lf) {
                while (ff <= lf) {
                    copyRecord(ff++, numRecordMerged++, b, sortedChunk);
                }
            } else {
                while(fs <= ls) {
                    copyRecord(fs++, numRecordMerged++, b, sortedChunk);
                }
            }
        }
        /*
        if(!checkSorted(sortedChunk)) {
            
            System.out.println("not sorted after merge, swich point is: " + getSwitchIndex(sortedChunk));
            //throw new RuntimeException("!!!!!!!! Not correctly sorted !!!!!!!!!");
        } else {
            System.out.println("The chunk is properly sorted after mergeSortedChunk");
        } 
        */
        raf.seek(pos);
        raf.write(sortedChunk);

        timer.stopTiming();
        System.out.println("Time used to merge a chunk: " + timer.getTimeUsed());
        return false;
    }
    private void copyRecord(int sourceIndex, int targetIndex, byte [] source, byte [] target) {
        sourceIndex *= RECORDSIZE;
        targetIndex *= RECORDSIZE;
        for(byte i = 0; i < RECORDSIZE; i++) {
            target[targetIndex++] = source[sourceIndex++];
        }
    } 
    /**
     * @return return true if the temp array is aready sorted before calling this function
     */
    private boolean sortChunk(RandomAccessFile raf, long pos, byte [] bytes) throws IOException {
        raf.seek(pos);
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
        raf.seek(pos);
        raf.write(bytes);
        return false;
    }
    /**
     * Heap sort algorithm is used
     * @param b - the byte array to be sorted
     * @param recordSize - the size of each record
     * @param startByte - the index of the first byte of the field to sort upon 
     * @param numBytes - number of bytes of the field to sort upon. e.g., 4 for int
     */
    private void sortRecords(byte [] bytes) {
        int numRecords = NUMBYTESPERCHUNK/RECORDSIZE;
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
        for(int i = 1; i < NUMRECORDSPERCHUNK; i++) {
            if(compare(i, i-1, bytes) < 0) {
                return false;
            }
        }
        return true;
    }
    // check if the sorted first half of the records in b are all
    // smaller or equal to the sorted second half of the records in b
    // return true if so, otherwise return false 
    private boolean isSorted(byte [] bytes) {
        // index of the first element of second half
        int firstSecondHalf = NUMRECORDSPERCHUNK/2; 
        // index of the last element of first half
        int lastFirstHalf = firstSecondHalf - 1;
        
        return (compare(lastFirstHalf, firstSecondHalf, bytes) != 1);
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
    private void outputPeptides() throws IOException {
        TimeUtils timer = new TimeUtils();
        timer.startTiming();
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream
                                     (new FileOutputStream(peptideIndexFile)));

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
        dos.close();
    }

    // index is the index for the fasta in the database
    private int outputPeptides(int index, Fasta f, DataOutputStream dos) throws IOException {
        
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
                    dos.writeInt(intMass); // mass
                    dos.writeInt(index);  // index of the fasta
                    dos.writeChar(i);  // start of the peptide
                    dos.writeChar(j); // end of the peptide
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
            ProcessDb se = new ProcessDb(fastaFile, sp);
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
}


