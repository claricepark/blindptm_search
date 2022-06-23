package edu.scripps.pms.util.io;

import edu.scripps.pms.util.seq.Fasta;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.HashSet;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.File;
//import edu.scripps.pms.util.dtaselect.Protein;
import java.util.ArrayList;
//import edu.scripps.pms.util.dtaselect.Peptide;


/**
 * @author  Tao Xu
 * @author  Robin Park
 * @version $Id: FastaReader.java,v 1.24 2014/07/08 23:52:13 rpark Exp $
 */
public class FastaReader
{

    public static final char FIRSTCHAROFDEFLINE = '>';
    public static final int DEFAULTSEQENCELENGTH = 200;
    public static final String RSTART = ">Reverse_"; // for reverse entries
    public static final String RRSTART = ">Reverse_Reverse_"; // for reverse_reverse entries
    public static final String USAGE = "USAGE: fastachecker fasta_file_or_folder [no_for_not_to_check_redundency]";
   

    // Becareful, might need lots of memory
    public static List <Fasta> getFastaList(InputStream is) throws IOException {
        List fastaList = new LinkedList();
        for (Iterator <Fasta> fastas = getFastas(is); fastas.hasNext();) {
            fastaList.add(fastas.next());
        }
        return fastaList;
    }
    public static Iterator<Fasta> getFastas(String fastaFileName) throws IOException {
        FileInputStream fis = new FileInputStream(fastaFileName);
        return getFastas(fis);
    } 
    public static Iterator <Fasta> getFastas(final InputStream is) throws IOException {
        return new Iterator() {

            private String lastLine = ""; // remember the last line read
            private BufferedReader br;

            {
                br = new BufferedReader(new InputStreamReader(is));
                // remove the potential empty lines and get the first defline
                while ((lastLine = br.readLine()) != null && lastLine.equals(""));

                if (lastLine.charAt(0) != FIRSTCHAROFDEFLINE) {
                    throw new FileFormatUnknownException();
                }
            }

            public boolean hasNext() {
                return lastLine != null;
            }

            public Object next() {

                Fasta fasta = null;
                try {
                    fasta = getFasta();
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return fasta;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private Fasta getFasta() throws IOException {

                StringBuffer sb = new StringBuffer(DEFAULTSEQENCELENGTH);
		String defline = lastLine;

                // if the line read is a empty string, ignore it
                while ( (lastLine = br.readLine()) != null && (lastLine.equals("")
                    || lastLine.charAt(0) != FIRSTCHAROFDEFLINE)) {
                    //System.out.println(lastLine);
                    if (!lastLine.equals("")) {
                        String line = lastLine.trim();
                        sb.append(line);
                    }
                }

                // the lastLine should be the defline
                // and sb.toString should be the sequence
                return new Fasta(defline, sb.toString());
            }

            protected void finalize() throws IOException {
                br.close();
                //System.out.println("Finalized");
            }
        };
    }
 
    public static int [] checkForwardReverseEntries(String fastafile) throws IOException {
        // for number of forward and reverse entries
        // first element for number of forward entries
        // second elements for number of entries starts with >Reverse_
        // third element for number of entries starts with >Reverse_Reverse_
        int [] numfr = new int [3]; 

        String line = ""; // remember the last line read
        BufferedReader br = null;
        try {

            br = new BufferedReader(new InputStreamReader(new FileInputStream(fastafile)));
            while (line != null) {
                if(line.length() > 0 && line.charAt(0) == '>') { // defline
                    if(line.startsWith(RSTART)) { // reverse entry
                        numfr[1]++; 
                    } else if(line.startsWith(RRSTART)) { // reverse_reverse entry
                        numfr[2]++;
                    } else { // regular entry
                        numfr[0]++;
                    }    
                } 
                line = br.readLine();
               
            }
        
        // remove the potential empty lines and get the first defline
        } catch(IOException e) { 
            e.printStackTrace();
        } finally {
            if (br != null) br.close();
        }
        return numfr;
        
    }

    public static void main(String args[]) throws IOException
    {
/*
        for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(args[0])); itr.hasNext(); ) { 
            Fasta fasta = itr.next();
            String defline = fasta.getDefline();
            if(defline.contains("Escherichia coli")) {
                System.out.println(">" + defline);
                System.out.println(fasta.getSequence());
            }
        }
*/
        if(args.length < 1) {
            System.out.println(USAGE);
            System.exit(1);
        }
         boolean checkRedundancy = args.length == 2? false : true;
         String dir = args[0];
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         if(currentDir.isDirectory()) {
            for(String s : currentDir.list()) {
                if(s.endsWith(".fasta") || s.endsWith(".FASTA")) {
                    files.add(dir + "/" + s);
                    //System.out.println(s);
                }
            }
         } else {
             files.add(dir);
         }

        for(int i = 0; i < files.size(); i++) {


            //String fastafile = args[0];
            String fastafile = files.get(i);
            int [] numfretnries = checkForwardReverseEntries(fastafile);
            //System.out.println("Number of forward entries: " + numfretnries[0]);
            //System.out.println("Number of reverse entries: " + numfretnries[1]);
            //System.out.println("Number of reverse_reverse entries: " + numfretnries[2]);
            String hasReverseReverse = numfretnries[2] > 0 ? " It contains Reverse_Reverse entries. " : "";
            System.out.println("");
            if(numfretnries[0] == numfretnries[1] && numfretnries[2] == 0) {
                System.out.println("Database " + fastafile + " IS properly reversed. " + hasReverseReverse + " # forward entriesis: " + numfretnries[0] + "\t# reversed entries: " + numfretnries[1] + "\t# Reverse_Reverse entries: " + numfretnries[2] );
            } else {
                System.out.println("Database " + fastafile + " IS NOT properly reversed.  # forward entriesis: " + numfretnries[0] + "\t# reversed entries: " + numfretnries[1] + "\t# Reverse_Reverse entries: " + numfretnries[2] );

                //System.out.println("Improperly reversed fasta database: " + fastafile);
            }
            if(checkRedundancy) {
                int numEntries = 0;
                HashSet<String> accessions = new HashSet<String>(1000000);
                HashSet<String> sequestLikeAccs = new HashSet<String>(1000000); 
                for (Iterator itr = FastaReader.getFastas(new FileInputStream(fastafile)); itr.hasNext(); ) {
                    Fasta fasta = (Fasta) itr.next();
        //System.out.println(fasta.getSequestLikeAccession());
                    numEntries++;
                
                    String defLine = fasta.getDefline();
                    String seq = fasta.getSequence();
                    String accession = fasta.getAccession(); 
                    if(accessions.contains(accession)) {
                        System.out.println("Duplicate accession: " + accession + "\tDefline: " + defLine);
                    }
                    accessions.add(accession);
            
                    String sequestlikeac = fasta.getSequestLikeAccession();
                    if(sequestlikeac.length() > 40) {
                        sequestlikeac = fasta.getSequestLikeAccession().substring(0, 41);
                    }
                    //System.out.println("acc: " + fasta.getAccession() + " seqacc: "  + fasta.getSequestLikeAccession() + " defline: " + fasta.getDefline()); 
                
                    if(sequestLikeAccs.contains(sequestlikeac)) {
                        System.out.println("Duplicate sequest like accession: " + sequestlikeac + "\tDefline: " + defLine);
                    }
                    sequestLikeAccs.add(sequestlikeac);
                }

                //System.out.println("In fasta file " + args[0] + ":");
                System.out.println("Number of protein entries: " + numEntries);
                System.out.println("Number of unique accessions: " + accessions.size());
                System.out.println("Number of unique SEQUEST like accessions: " + sequestLikeAccs.size());
        
            
        /*
                for (Iterator itr = FastaReader.getFastas(new FileInputStream(args[0])); itr.hasNext(); ) {
                    Fasta fasta = (Fasta) itr.next();
                    String defLine = fasta.getDefline();
                    String seq = fasta.getSequence();
                    if(defLine.startsWith("Rever")) {
                        //System.out.println("Reversed: " + defLine);
                        if(seq.endsWith("M")) {
                            seq = seq.substring(0, seq.length() -1);
                        } 
                    } else {
                        if(seq.startsWith("M")) {
                            seq = seq.substring(1, seq.length());

                        }
                        //System.out.println("Regular: " + defLine);
                    }
                    
                System.out.println(">" + defLine);
                System.out.println(seq);
                }
        */      

            }
        }
    }
}
