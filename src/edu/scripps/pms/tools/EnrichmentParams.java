package edu.scripps.pms.tools;

import java.io.*;
import java.util.*;

/* EnrichmentParams
   Initiated 080105
   John Venable, The Scripps Research Institute

   * Reads in a .params file for enrichmentcalc settings
   */

class EnrichmentParams {
	public int  resolution = 100000;
	public int  shift = 5;
	public int  enrichmax = 100;
	public int  enrichmin = 80;
	public float  enrichdelta = 0.1f;

	public EnrichmentParams() {
                this.readParams();
        }

	public void readParams() {
                try {
                        //defines local variables and checks to see that the precursorcalc.params file exists
                        File            ParamsFile = new File("enrichmentcalc.params");
                        if (ParamsFile.exists()) {
                                FileReader      InputFileReader = new FileReader(ParamsFile);
                                BufferedReader  Incoming = new BufferedReader(InputFileReader);
                                String          LineBuffer;
                                String          WholeLine;
                                StringTokenizer Parser;
                                WholeLine = Incoming.readLine();
                                //Reads in each line in enrichmentcalc.params and assigns values to appropriate instance variables
                                while (WholeLine != null) {
                                        Parser = new StringTokenizer(WholeLine, "\t");
                                        if (Parser.hasMoreTokens()) {
                                                LineBuffer = Parser.nextToken();
                                                if (LineBuffer.startsWith("#")) {
                                                        // It's a comment; ignore it.
                                                }
                                                else if (LineBuffer.equals("Resolution")) {
                                                        this.resolution = new Integer(Parser.nextToken()).intValue();
                                                }
                                                else if (LineBuffer.equals("Shifts_Allowed")) {
                                                        this.shift = new Integer(Parser.nextToken()).intValue();
                                                }
						else if (LineBuffer.equals("Enrichment_Max")) {
                                                        this.enrichmax = new Integer(Parser.nextToken()).intValue();
                                                }
						else if (LineBuffer.equals("Enrichment_Min")) {
                                                        this.enrichmin = new Integer(Parser.nextToken()).intValue();
                                                }
						else if (LineBuffer.equals("Enrichment_Delta")) {
                                                        this.enrichdelta = new Float(Parser.nextToken()).floatValue();
                                                }
                   				else {
                                                        System.out.println("I don't understand this option in enrichment.params.");
                                                        System.out.println(WholeLine);
                                                        System.exit(0);
                                                }
                                        }
                                        WholeLine = Incoming.readLine();
                                }
                                // Close file
                                Incoming.close();
                        } else {
				System.out.println("enrichmentcalc.params file not found");
			}
                }
                //error handling
                catch (IOException failure) {
                        System.out.println("Error while reading enrichmentcalc.params");
                        System.out.println(failure.toString());
                }
                //return this;
        }

        public void printParams() {
                System.out.println("EnrichmentCalc Parameters");
                System.out.println("Resolution " + this.resolution);
                System.out.println("Shifts Allowed " + this.shift);
		System.out.println("Enrichment_Max " + this.enrichmax);
                System.out.println("Enrichment_Min " + this.enrichmin);
		System.out.println("Enrichment_Delta " + this.enrichdelta);
        }
}
