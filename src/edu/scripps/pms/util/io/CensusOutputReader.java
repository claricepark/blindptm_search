/**
 * @file CensusOutputReader.java
 * This is the source file for CensusOutputReader
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.util.io;

import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.census.CensusPeptide; 
import edu.scripps.pms.util.census.CensusProtein; 

import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for calculating the z-score for delta mass 
public class CensusOutputReader {
    public static final String USAGE = "\n\n!!! USAGE: java CensusOutputReader -i input -d downthreshold -u upthreshold!!!";
    private static String inputFile = null;
    private static double upthreshold = 1.5;
    private static double downthreshold = 0.75;
    private static boolean modifiedOnly = true;

    private String file;
    private ArrayList<CensusProtein> proteins = new ArrayList<CensusProtein>(10000);

    public CensusOutputReader(String file) throws IOException {
        this.file = file;
        readProteins();  
    }

    public Iterator<CensusProtein> getCensusProteins() {
        return proteins.iterator();
    }
    private void readProteins() throws IOException {
        BufferedReader  br = new BufferedReader(new FileReader(file));
        String lastLine = br.readLine();

        while ((lastLine = br.readLine()) != null && (lastLine.startsWith("H"))); 
        while ((!lastLine.startsWith("P")) && (lastLine = br.readLine()) != null); 
        if(!lastLine.startsWith("P")) {
            System.err.println("The census output file format is probably wrong, check if the first PLine is properly formated");
            System.exit(1);

        }

        while(lastLine != null) {
            CensusProtein p = new CensusProtein(lastLine);
            boolean isnewprotein = true;
            proteins.add(p);
            while ((lastLine = br.readLine()) != null) {
                if(lastLine.startsWith("P")) {
                    if(!isnewprotein) {
                        break;
                    } 
                    p.addProteinLine(lastLine); 
                }else if(lastLine.startsWith("S")) {
                    p.addPeptide(lastLine);
                    isnewprotein = false;
                } else {
                    // ignore unexpected lines
                }
            }
            
        } 
        br.close();
    }
    public static void main(String args[]) throws Exception {

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - census output file");
        Option upOpt = new Option
            ("u", "upthreshold", true,  "threshold for upregulated, e.g. 1.5");
        Option downOpt = new Option
            ("d", "downthreshold", true,  "threshold for downregulated,e.g., 0.75");
        //inputOpt.setRequired(true);
        //upOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(upOpt);
        opts.addOption(downOpt);

        BasicParser cliParser = new BasicParser();
        CommandLine cli = null;

        try {
            cli = cliParser.parse(opts, args);
            
        } catch (Exception e) {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp(USAGE, opts);
            e.printStackTrace();
            System.out.println
                ("\nPlease use \"java -Xmx1500M ... \" when processing a large dataset.");
            System.exit(1);
        }
        inputFile = cli.getOptionValue("i");

        if(cli.getOptionValue("u") != null) {
            upthreshold = Double.parseDouble(cli.getOptionValue("u"));        
        }
        if(cli.getOptionValue("d") != null) {
            downthreshold = Double.parseDouble(cli.getOptionValue("d"));        
        }

        try {
            CensusOutputReader cor = new CensusOutputReader(inputFile);
            for(Iterator<CensusProtein> it = cor.getCensusProteins(); it.hasNext();) {
                CensusProtein p = it.next();
                System.out.println(p.getRepresentativeAccession() + "\t" + p.getProteinRatio() + "\t" + p.getMaxRatioPeptide().getRatio() + "\t" + p.getMinRatioPeptide().getRatio() +  "\t" + p.getDescription());
            }
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }
}
