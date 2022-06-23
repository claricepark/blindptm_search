package edu.scripps.pms.tools;

import java.io.BufferedReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.PeakList;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class SpectrumValidator
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide = true;
    private int lineNumber = 1;


    public static void main(String args[]) throws IOException
    {
        String eachFile;

        for (Iterator fileItr = getAllFiles(args[0]); fileItr.hasNext(); ) {
            eachFile = fileItr.next().toString();

            try {
                if (eachFile.endsWith("ms2")) {
                    System.out.println("Validating file " + eachFile + "...");

                    String line = null;
                    SpectrumReader sr = new SpectrumReader(args[0] + File.separator + eachFile, "ms2");
                    Iterator<PeakList> it = sr.getSpectra();
                    int counter = 0;
                    int numPeaks = 0;
                    boolean sortByIntensity = true;

                    for (Iterator itr = sr.getSpectra(); itr.hasNext(); ) {
                        PeakList pkList = (PeakList) itr.next();

                        pkList.getListType();
                        pkList.getLoscan();
                        pkList.getHiscan();
                    }

                    System.out.println("");
                }
            } catch (Exception e) {
                System.out.println("error : " + e);
            }
        }

        System.out.println("Validating sqt files finished");
    }


    public static Iterator getAllFiles(String path) {
        File f = new File(path);
        return Arrays.asList(f.list()).iterator();
    }


}

