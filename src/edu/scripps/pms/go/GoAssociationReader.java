/**
 * @file GoAssociationReader.java
 * This is the source file for edu.scripps.pms.util.spectrum.GoAssociationReader
 * @author Tao Xu
 * @date $Date: 2006/07/13 00:53:43 $
 */



package edu.scripps.pms.go;

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.HashSet;
import java.util.Collections;

public class GoAssociationReader {


    // delimiter of m2z and intensity
    // the format of ms data file, e.g., ms1, ms2, dta

    // if the spectraDelimiter is "S", then isNewFormat is true,
    // if the spectraDelimiter is ":", then isNewFormat is false
    // delimiter for spectra, could be ':' or 'S', depends on the ms file type
    private String associationFile;
    public GoAssociationReader (String fileName) throws IOException {
        associationFile = fileName;
    }

    public static void main(String args[]) throws Exception {
        HashSet<String> ipiAccs = new HashSet<String>();
        HashSet<String> swissProtAccs = new HashSet<String>();
        HashSet<String> taxons = new HashSet<String>();
        GoAssociationReader gr = new GoAssociationReader(args[0]);

        Iterator<GoAssociation> it = gr.getGoAssociations();
        int counter = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            GoAssociation g = it.next();
            swissProtAccs.add(g.getDbObjectId());
            ipiAccs.add(g.getDbObjectSynonym());
            taxons.add(g.getTaxon());
            counter++;
            //if (counter > 20)
            //break;
//            System.out.println("IPI: " + g.getDbObjectSynonym() + "\tSwissProt: " + g.getDbObjectId());
        }
        System.out.println("Total number of associations processed: " + counter);
        System.out.println("numIPI: " + ipiAccs.size() + "\tnumSwissProt: " + swissProtAccs.size() + "\tnumTaxon: " + taxons.size());

        System.out.println("Finished");
    }


    public Iterator <GoAssociation> getGoAssociations() throws IOException {
        
        final BufferedReader br = new BufferedReader(new FileReader(associationFile));
        return new Iterator() {
            String lastLine = null;
            public boolean hasNext() {
                boolean hasNext = true;
                try {
                    lastLine = br.readLine();
                    while(lastLine != null && lastLine.startsWith("!")) {
                        lastLine = br.readLine();
                    } 
                    hasNext = (lastLine != null);
                    if(!hasNext) {
                        br.close(); 
                    } 
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return hasNext;
            }

            public Object next() {
                return new GoAssociation(lastLine);
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }
            private void closeDataFile() throws IOException {
                br.close();
            }
        };
    }




}
