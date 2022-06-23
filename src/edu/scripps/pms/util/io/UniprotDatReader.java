package edu.scripps.pms.util.io;

import edu.scripps.pms.util.seq.UniProtProtein;
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
import edu.scripps.pms.util.dtaselect.Protein;
import java.util.ArrayList;
import edu.scripps.pms.util.dtaselect.Peptide;

public class UniprotDatReader
{

    public static final int DEFAULTSEQENCELENGTH = 1000;

    public static Iterator<UniProtProtein> getUniProtProteins(String fastaFileName) throws IOException {
        FileInputStream fis = new FileInputStream(fastaFileName);
        return getUniProtProteins(fis);
    } 
    public static Iterator <UniProtProtein> getUniProtProteins(final InputStream is) throws IOException {
        return new Iterator() {

            private String lastLine = ""; // remember the last line read
            private BufferedReader br;

            {
                br = new BufferedReader(new InputStreamReader(is));
                // remove the potential empty lines and get the first defline
                while ((lastLine = br.readLine()) != null && lastLine.equals(""));

                if (!lastLine.startsWith("ID")) {
                    throw new FileFormatUnknownException();
                }
            }

            public boolean hasNext() {
                return lastLine != null;
            }

            public Object next() {

                UniProtProtein fasta = null;
                try {
                    fasta = getUniProtProtein();
                } catch (IOException e) {
                    e.printStackTrace();
                }

                return fasta;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private UniProtProtein getUniProtProtein() throws IOException {

                StringBuffer sb = new StringBuffer(DEFAULTSEQENCELENGTH);
		String defline = lastLine;
//System.out.println(lastLine);
                UniProtProtein protein = new UniProtProtein(lastLine);
                // if the line read is a empty string, ignore it
                while ( (lastLine = br.readLine()) != null && (lastLine.equals("")
                    || !lastLine.startsWith("ID"))) {
                    //System.out.println(lastLine);
                    if (!lastLine.equals("")) {
                        //String line = lastLine.trim();
                        //sb.append(line);
                        protein.addLine(lastLine);
                    }
                }

                // the lastLine should be the defline
                // and sb.toString should be the sequence
                
                return protein; 
            }

            protected void finalize() throws IOException {
                br.close();
                //System.out.println("Finalized");
            }
        };
    }


    public static void main(String args[]) throws IOException {
        String uniprotdatfile = "/data/6/xmhan/database/uniprot_sprot_human_11-14-2007.dat";
        String salivaUniProtProtein = "/data/6/xmhan/database/EBI-IPI_human_3.37_12-05-2007_original.fasta";
        
        for (Iterator<UniProtProtein> itr = UniprotDatReader.getUniProtProteins(new FileInputStream(uniprotdatfile)); itr.hasNext(); ) { 
            UniProtProtein fasta = itr.next();
            
        //    String defline = fasta.getDefline();
            //if(defline.contains("Escherichia coli")) {
          //      System.out.println(">" + defline);
           //     System.out.println(fasta.getSequence());
            //}
        }
    }
}
