package edu.scripps.pms.qcorr;

import java.io.*;
import java.util.*;

/* ParamsFile
   Initiated 080105
   John Venable, The Scripps Research Institute

   * Reads in a .params file for QCorr settings
   */

class ParamsFile {
        //defines instance variables
	
        public int  NumIons = 20;
        public float  InjWindow = 3f;
        public float MassAccuracy = 0.1f;
        public boolean  Normalize = false;
	public boolean  FilterPeaks = false;
	public int  PeakThreshold = 200;
	public boolean FilterSpectra = false;
	public boolean  OverWriteMS2 = false;
	public boolean  ReplaceMW = false;
	public boolean  AssignCharge = false;
	public float[] FilterThreshold = new float[4];	
        public boolean  OutputFile = false;
        public boolean OutputCCArray = false;
	public boolean DeChargeSpectrum = false;
	public float IsotopeAccuracy = 500;//units of ppm
	public boolean AssignMWFromMS1 = false;
	private String homedir;


        public ParamsFile(String homedir) {
        	this.homedir = homedir;
		this.readParams(homedir);
    	}

	//initialize routine reads in "precursorcalc.params" file and sets variables to appropriate values
        public void readParams(String homedir) {
                try {
                        //defines local variables and checks to see that the precursorcalc.params file exists
			System.out.println("Params Directory = " + homedir);
                        File            ParamsFile = new File("QCorr.params");
			File		HomeParamsFile = new File(homedir + "QCorr.params"); 		
                        FileReader      InputFileReader = null;
			if (ParamsFile.exists()) {
                                InputFileReader = new FileReader(ParamsFile);
			} else if (HomeParamsFile.exists()) {
				InputFileReader = new FileReader(HomeParamsFile);
			}
			if (ParamsFile.exists() || HomeParamsFile.exists()) {	
				BufferedReader  Incoming = new BufferedReader(InputFileReader);
                                String          LineBuffer;
                                String          WholeLine;
				String		temp;
				float		p1;
				float           p2;
				float           p3;
				float		p4;
                                StringTokenizer Parser;
				StringTokenizer Threshold;
                                WholeLine = Incoming.readLine();
                                //Reads in each line in precursorcalc.params and assigns values to appropriate instance variables
                                while (WholeLine != null) {
                                        Parser = new StringTokenizer(WholeLine, "\t");
                                        if (Parser.hasMoreTokens()) {
                                                LineBuffer = Parser.nextToken();
                                                if (LineBuffer.startsWith("#")) {
                                                        // It's a comment; ignore it.
                                                }
						else if (LineBuffer.equals("NumIons")) {
                                                        this.NumIons = new Integer(Parser.nextToken()).intValue();
                                                }
                                                else if (LineBuffer.equals("InjWindow")) {
                                                        this.InjWindow = new Float(Parser.nextToken()).floatValue();
                                                }
						else if (LineBuffer.equals("MassAccuracy")) {
                                                        this.MassAccuracy = new Float(Parser.nextToken()).floatValue();
                                                }
						else if (LineBuffer.equals("Normalize")) {
                                                        this.Normalize = new Boolean(Parser.nextToken()).booleanValue();
                                                }	
                                                else if (LineBuffer.equals("FilterPeaks")) {
                                                        this.FilterPeaks = new Boolean(Parser.nextToken()).booleanValue();
                                                }
                                                else if (LineBuffer.equals("PeakThreshold")) {
                                                        this.PeakThreshold = new Integer(Parser.nextToken()).intValue();
                                                }
                                                else if (LineBuffer.equals("FilterSpectra")) {
                                                        this.FilterSpectra = new Boolean(Parser.nextToken()).booleanValue();
                                                }
                                                else if (LineBuffer.equals("OverWriteMS2")) {
                                                        this.OverWriteMS2 = new Boolean(Parser.nextToken()).booleanValue();
                                                }
						else if (LineBuffer.equals("ReplaceMW")) {
                                                        this.ReplaceMW = new Boolean(Parser.nextToken()).booleanValue();
                                                }
						else if (LineBuffer.equals("AssignCharge")) {
                                                        this.AssignCharge = new Boolean(Parser.nextToken()).booleanValue();
                                                }
                                                else if (LineBuffer.equals("FilterThreshold")) {
							temp = Parser.nextToken();
							Threshold = new StringTokenizer(temp, ",");
							p1 = new Float(Threshold.nextToken()).floatValue();
							p2 = new Float(Threshold.nextToken()).floatValue();
							p3 = new Float(Threshold.nextToken()).floatValue();
							p4 = new Float(Threshold.nextToken()).floatValue();
                                                        this.FilterThreshold[0] = p1;
							this.FilterThreshold[1] = p2;
							this.FilterThreshold[2] = p3;
							this.FilterThreshold[2] = p4;
                                                }
						else if (LineBuffer.equals("OutputFile")) {
                                                        this.OutputFile = new Boolean(Parser.nextToken()).booleanValue();
                                                }
                                                else if (LineBuffer.equals("OutputCCArray")) {
                                                        this.OutputCCArray = new Boolean(Parser.nextToken()).booleanValue();
                                                }
						else if (LineBuffer.equals("DeChargeSpectrum")) {
                                                        this.DeChargeSpectrum = new Boolean(Parser.nextToken()).booleanValue();
                                                }
						else if (LineBuffer.equals("IsotopeAccuracy")) {
                                                        this.IsotopeAccuracy = new Float(Parser.nextToken()).floatValue();
                                                }
						else if (LineBuffer.equals("AssignMWFromMS1")) {
                                                        this.AssignMWFromMS1 = new Boolean(Parser.nextToken()).booleanValue();
                                                }
                                                else {
                                                        System.out.println("I don't understand this option in QCorr.params.");
                                                        System.out.println(WholeLine);
                                                        System.exit(0);
                                                }
                                        }
                                        WholeLine = Incoming.readLine();
                                }
                                // Close file
                                Incoming.close();
                        } else {
				System.out.println("QCorr.params does not exist! Using default criteria.");
			}
                }
                //error handling
                catch (IOException failure) {
                        System.out.println("Error while reading precursorcalc.params");
                        System.out.println(failure.toString());
                }
                //return this;
        }

	public void printParams() {
		System.out.println("QCorr Parameters");
		System.out.println("NumIons " + this.NumIons);
        	System.out.println("InjWindow " + this.InjWindow);
        	System.out.println("MassAccuracy " + this.MassAccuracy);
        	System.out.println("Normalize " + this.Normalize);
        	System.out.println("Filterpeaks " + this.FilterPeaks);
        	System.out.println("PeakThreshold " + this.PeakThreshold);
        	System.out.println("FilterSpectra " + this.FilterSpectra);
        	System.out.println("OverWriteMS2 " + this.OverWriteMS2);
		System.out.println("ReplaceMW " + this.ReplaceMW);
        	System.out.println("AssignCharge " + this.AssignCharge);
        	System.out.println("FilterThreshold1 " + this.FilterThreshold[0]);
		System.out.println("FilterThreshold2 " + this.FilterThreshold[1]);
		System.out.println("FilterThreshold3 " + this.FilterThreshold[2]);
		System.out.println("FilterThreshold4 " + this.FilterThreshold[3]);
        	System.out.println("OutputFile " + this.OutputFile);
        	System.out.println("OutputCCArray " + this.OutputCCArray);
		System.out.println("AssignMWFromMS1 " + this.AssignMWFromMS1);
	}
}

