package edu.scripps.pms.msfilter;

import java.io.*;
import java.util.*;
import java.text.*;
import java.lang.*;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.SpectrumIndexReader;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;



/*
    Test program for fragment ion intensity profile.
    Daniel Cociorva and John Venable, September 12, 2004.
*/

public class MSFilter {
    public String             MS2Filename = null;
    public String	      NewMS2Filename = null;
    public IniFile            Config = new IniFile();
    public int		      LastSpectrum = 0;
    public Spectrum           SpectrumList = new Spectrum();
    public SpectralMatrix     Matrix = new SpectralMatrix();

    public static void main(String args[]) throws Exception {
        MSFilter            AppObject = new MSFilter();
	int argcntr = 0;

	if (args.length > 0) {
		AppObject.Go(args);
	} else {
		System.out.println("Must enter an ms2 file...");
	}
    }

    public void Go(String[] args) throws Exception {
	File            currentDirectory = new File(System.getProperty("user.dir"));
	EndsWithFilter  uniMS2Filter = new EndsWithFilter(".ms2");
	String[]        uniMS2List = currentDirectory.list(uniMS2Filter);
	int             dirIterator;
	String          currentMS2File;
	PeakList        list;
	String[]        data;
	ArrayList<String>       spectrumData;
	int             specCount = 1;
	int             minPeaks = 10;
	Runtime         systemCommand = Runtime.getRuntime();
	Process         process;
	List<String>    hLines = new ArrayList();
        int 		SpecCount = 0;
	int		middleSpectrum = 0;
	boolean		singleSpectrum = false;
	boolean 	quit = false;

	//Read in params file
	Config.Initialize();
	
	if (args.length == 2 && new Integer(args[1]).intValue() > 0) {
		middleSpectrum = new Integer(args[1]).intValue();
		singleSpectrum = true;
		currentMS2File = args[0];
		System.out.println("Reading " + currentMS2File + "...");
		SpectrumReader ms2 = new SpectrumReader(currentMS2File, "ms2");
		Iterator<PeakList> spectrumList = ms2.getSpectra();
		Spectrum spectrum = new Spectrum();
		while (spectrumList.hasNext() && !quit) {
			list = spectrumList.next();
			System.out.println("Looking for spectrum " + middleSpectrum + "...");
			data = spectrum.processSpectrum(singleSpectrum, Config, list, specCount, currentMS2File, middleSpectrum);
			if(data[0] == "quit") {
				quit = true;
			}
			specCount++;
		}

	} else {
		//read in MS2 and Parse through MS2. Each spectrum is processed in turn
		for(dirIterator = 0; dirIterator < args.length; dirIterator++) {
			currentMS2File = args[dirIterator];
			System.out.println("Reading " + currentMS2File + "...");
			SpectrumReader ms2 = new SpectrumReader(currentMS2File, "ms2");
			Iterator<PeakList> spectrumList = ms2.getSpectra();
			Spectrum spectrum = new Spectrum();
			while (spectrumList.hasNext()) {
				list = spectrumList.next();
				System.out.println("Processing spectrum " + list.getLoscan() + "...");
				data = spectrum.processSpectrum(singleSpectrum, Config, list, specCount, currentMS2File, middleSpectrum);
				specCount++;
			}
		}
	}
    }
}

class IniFile {
        public float   FragmentTolerance = 0f;
	public int     NumScans = 0;
	public int     ScanInterval = 0;
	public float   Threshold = 0f;
	public int     MaxPeaks = 200;
	public void    Initialize() {
            // Set up a list of residues
            try {
                File            MSFilterIni = new File("MSFilter.ini");
                if (MSFilterIni.exists()) {
                    FileReader      InputFileReader = new FileReader(MSFilterIni);
                    BufferedReader  Incoming = new BufferedReader(InputFileReader);
                    String          LineBuffer;
                    String          WholeLine;
                    StringTokenizer Parser;
                    WholeLine = Incoming.readLine();
                    while (WholeLine != null) {
                        Parser = new StringTokenizer(WholeLine);
                        if (Parser.hasMoreTokens()) {
                            LineBuffer = Parser.nextToken();
                            if (LineBuffer.startsWith("#")) {
                                // It's a comment; ignore it.
                            }
			    else if (LineBuffer.equals("MaxPeaks")) {
                                this.MaxPeaks = new Integer(Parser.nextToken()).intValue();
                            }
			    else if (LineBuffer.equals("Correlation")) {
                                this.Threshold = new Float(Parser.nextToken()).floatValue();
                            }
                            else if (LineBuffer.equals("FragmentTolerance")) {
                                this.FragmentTolerance = new Float(Parser.nextToken()).floatValue();
                            }
			    else if (LineBuffer.equals("NumberScans")) {
                                this.NumScans = new Integer(Parser.nextToken()).intValue();
                            }
			    else if (LineBuffer.equals("ScanInterval")) {
                                this.ScanInterval = new Integer(Parser.nextToken()).intValue();
                            }
                            else {
                                System.out.println("I don't understand this option in MSFilter.ini:");
                                System.out.println(WholeLine);
                                System.exit(0);
                            }
                        }
                        WholeLine = Incoming.readLine();
                    }
                    // Close file
                    Incoming.close();
                } else {
			System.out.println("Error while reading MSFilter.ini");
			System.exit(0);
		}
            }
            catch (IOException failure) {
                System.out.println("Error while reading MSFilter.ini");
                System.out.println(failure.toString());
            }
        }
}

class SpectrumNumbers {
	public int ScanNum;
	public SpectrumNumbers Next;
	public SpectrumNumbers Prev;
}

class Spectrum {
	public String[]      Charge;
	public String	     ZLines;
	public String	     SLine;
	public String	     DLine;
	public String	     ILine;
 	public double[]       MPlusH;
        public int           NPeaks = 0;
	public String        SpecNum;
	public double	     PrecursorMZ = 0f;
	public int	     NumberBY = 0;
	public double[]	     BandY_MZ;
	public String	     Sequence = "";
        public Spectrum      Prev;
        public Spectrum      Next;
	//public Isotoper IsotopeFilter = new Isotoper();
	public int           PeakCount = 0;
        public int           PostPeakCount = 0;
	public int      NumChargeStates = 0;
	public double	MaxInt = 0f;
	SpectralMatrix	Matrix = new SpectralMatrix();
	LinkedList<PeakList>	spectrumList = new LinkedList();
	ArrayList<Integer>	spectrumArray = new ArrayList();
	int		tCntr = 0;

	String[] processSpectrum(boolean singleSpectrum, IniFile params, PeakList list, int specCount, String currentMS2File, int middleSpectrum) {
		String[]	data;
		int		totalSpectra;
		PeakList	tList;
		int		iiloop;
		ListIterator<PeakList>	tempList;
		boolean		collectingData = false;
        	int           SpecCount = 0;
		int 		nclusters = 5; 
		int 		nrows;
		int 		ncolumns; 
		int 		npass = 100; 
		int[] 		clusterid;
		int		clusterScore;
		int		iilooper;
		float		preMZ;
		
		data = new String[2];
		totalSpectra = params.ScanInterval * params.NumScans;

		//add PeakList to LinkedList
		spectrumList.add(list);
		preMZ = (float)list.getPrecursorMass();
		spectrumArray.add(new Integer(list.getLoscan()));		

		if (!singleSpectrum && !collectingData) {
			collectingData = true;
			middleSpectrum = spectrumArray.get(tCntr).intValue();
			System.out.println("MiddleSpectrum = " + middleSpectrum + "\t" + tCntr);
			//if (list.getLoscan() == middleSpectrum) {
				System.out.println("Initializing Single Spectrum: " + middleSpectrum);
                                Matrix.initializeMatrix(params.NumScans*2+1, 0, list);
			//}
		}
		if (middleSpectrum > 0 && singleSpectrum) {
			if (list.getLoscan() == middleSpectrum) {
				System.out.println("Initializing Single Spectrum: " + middleSpectrum);
				Matrix.initializeMatrix(params.NumScans*2+1, 0, list);
			}
		}
		if (list.getLoscan() == middleSpectrum + totalSpectra) {
			collectingData = false;
                	System.out.println("Generating Spectrum Matrix For Scan: " + middleSpectrum);
			//load spectra before middle spectrum
                	for (iiloop=middleSpectrum-totalSpectra;iiloop<=middleSpectrum-params.ScanInterval;iiloop += params.ScanInterval) {
				tempList = spectrumList.listIterator(0);
                        	while (tempList.hasNext()) {
					tList = tempList.next();
                                	if (tList.getLoscan() == iiloop) {
                                        	SpecCount++;
                                                System.out.println("Loading Scans Before Midle Spectrum : " + tList.getLoscan() + "\t" + tList.numPeaks() + "\t" + SpecCount);
                                        	Matrix.LoadSpectrum(params, tList, SpecCount);
                                	}
                        	}
                	}

			//skip over middle spectrum
			//SpecCount++;
                	//load spectra after middle spectrum
			for (iiloop=middleSpectrum+params.ScanInterval;iiloop<=middleSpectrum+totalSpectra;iiloop += params.ScanInterval) {
				tempList = spectrumList.listIterator(0);
                        	while (tempList.hasNext()) {
                                	tList = tempList.next();
                                	if (tList.getLoscan() == iiloop) {
                                        	SpecCount++;
                                        	System.out.println("Loading Scans After Middle Spectrum: " + tList.getLoscan() + "\t" + tList.numPeaks() + "\t" + SpecCount);
                                        	Matrix.LoadSpectrum(params, tList, SpecCount);
                                        }
                        	}
                	}
			
			//Daniels code goes here//
			nrows = Matrix.NPeaks;
			clusterid = new int[nrows];
			ncolumns = params.NumScans*2+1;
			//cluster peaks
			//clusterScore = Matrix.clusterPeaks(nclusters, nrows, ncolumns, npass, clusterid);
			Matrix.Correlation(params.NumScans*2+1);
			//if (MakeProfiles){
                              //  Matrix.WriteOutputLine("Profiles_" + MS2Filename, SRunner.SLine, true);
                              // 	Matrix.WriteOutputLine("Profiles_" + MS2Filename, SRunner.ZLines, true);
                        	Matrix.WriteProfiles(params.NumScans, "Profiles_" + currentMS2File, true);
				Matrix.WriteOutputLine("Subset_" + currentMS2File, "S\t" + middleSpectrum + "\t" + middleSpectrum + "\t" + preMZ  + "\r\n", true);
				Matrix.WriteOutputLine("Remaining_" + currentMS2File, "S\t" + middleSpectrum + "\t" + middleSpectrum + "\t" + preMZ + "\r\n", true);
				int numZlines = list.getNumZlines();
				Iterator<Zline>	zLines = list.getZlines();
				while (zLines.hasNext()) {
					Zline	tzLine = zLines.next();
					Matrix.WriteOutputLine("Subset_"+ currentMS2File, "Z\t" + tzLine.getChargeState() + "\t" + tzLine.getM2z() + "\r\n", true);
					Matrix.WriteOutputLine("Remaining_"+ currentMS2File, "Z\t" + tzLine.getChargeState() + "\t" + tzLine.getM2z() + "\r\n", true);
				}
				Matrix.WriteSpectrum(0, "Remaining_" + currentMS2File, true, params, 0);
				Matrix.WriteSpectrum(0, "Subset_"+ currentMS2File, true, params, 1);

                        //}	
			//Matrix.printProfiles(SpecCount, clusterid);
                	System.out.println("Correlation Finished For Scan :\t" + middleSpectrum);
			tCntr++;
                }
		
		//Remove PeakLists from Linked List if their too far away from middleSpectrum
		if(spectrumList.getFirst().getLoscan() < middleSpectrum - totalSpectra) { 
			//System.out.println("Removing Scan :" + spectrumList.getFirst().getLoscan());
			tList = spectrumList.removeFirst();
		}
		
		//quit if only looking for single spectrum
		if(singleSpectrum && list.getLoscan() > middleSpectrum + totalSpectra) {
			data[0] = "quit";
		}
		return data;
	}

}
    
class SpectralMatrix {

        public int           NSpectra = 0;
        public int           NPeaks = 0;
        public int           SpecIndex = 0;
        public int           TopPeakIndex = 0;
	public double	     TrueAverage = 0f;
	public double         FalseAverage = 0f;
	public int	     TrueNum = 0;
	public int	     FalseNum = 0;
        public double[]       MZValues;
        public double[][]       RValues;
        public int[]         NPositive;
        public double[][]     Intensities;
	public boolean[]     ValidatedPeak;
	boolean correlateAgainstMaster = true;
        boolean correlateAgainstAll = false;
	public float		targetMZ = 761.7195f;

	public double	     Maxi = 0f;

	public void initializeMatrix(int n, int specIndex, PeakList list) {
		Peak	tPeak; 
		int	pIndex = 0;
		int 	nPeaks;
		
		nPeaks = list.numPeaks();		
		NPeaks = list.numPeaks();
		MZValues = new double[nPeaks];
            	RValues = new double[nPeaks][nPeaks];
		Intensities = new double[nPeaks][n];
		Iterator<Peak> peaks = list.getPeaks();	
		while(peaks.hasNext()) {
                	tPeak = peaks.next();
			MZValues[pIndex] = tPeak.getM2z();
			Intensities[pIndex][specIndex] = tPeak.getIntensity();
			if (correlateAgainstMaster && !correlateAgainstAll) {
				//determines peak to correlate against all others
			/*	if (list.get == 1) {
                        		if (PRunner.Intensity > Maxi) {
                                		Maxi = PRunner.Intensity;
                                		TopPeakIndex = Pindex;
                        		}
                		} else {
                        		if (PRunner.Intensity > Maxi && PRunner.MZ > MySpec.PrecursorMZ) {
                                		Maxi = PRunner.Intensity;
                                		TopPeakIndex = Pindex;
                        		}
                		}*/
				System.out.println(tPeak.getM2z() + "\t" + targetMZ);
				if (tPeak.getM2z() == targetMZ) {
					System.out.println("Found " + targetMZ);
					TopPeakIndex = pIndex;
				}
			}
			pIndex++;
		}
	}

	public int clusterPeaks(int nclusters, int nrows, int ncolumns, int npass, int clusterid[]) {
		int 	clusterScore = 0;
		//System.out.println("clusters " + nclusters+"\t"+nrows+"\t"+ncolumns+"\t"+clusterid.length+"\t"+Intensities[0][0]);
		clusterScore =Cluster.kcluster(nclusters, nrows, ncolumns, Intensities, npass, clusterid);
		return clusterScore;
	}

        public double TwoPeakCorrelation (int Pindex_x, int Pindex_y, int NSpectra) {
            double r = 0;
            int Sindex;
            double Sumx, Sumy, Sumxy, Sumx2, Sumy2, x, y;

            Sumx = 0; Sumy = 0; Sumx2 = 0; Sumy2 = 0; Sumxy = 0;
            for (Sindex = 0; Sindex < NSpectra; Sindex++) {
                x = Intensities[Pindex_x][Sindex];
                y = Intensities[Pindex_y][Sindex];
                Sumx += x;
                Sumy += y;
                Sumx2 += x*x;
                Sumy2 += y*y;
                Sumxy += x*y;
            }
            r = (Sumxy - Sumx*Sumy/NSpectra)/((double) Math.sqrt((Sumx2 - Sumx*Sumx/NSpectra) * (Sumy2 - Sumy*Sumy/NSpectra)));
            return r;
        }

	public void printProfiles(int nSpectra, int[] clusterid) {
		int Sindex;
		int Pindex;
		String Info;
		System.out.println("Starting to print profiles");
		//print scan header
		Info = "";
		for (Sindex = 0;Sindex < nSpectra; Sindex++) {
                        Info += ("\t" + Sindex);
                }
		System.out.println(Info);
		for (Pindex = 0; Pindex < NPeaks; Pindex++) {
			Info = "" + MZValues[Pindex];
			for (Sindex = 0;Sindex < nSpectra; Sindex++) {
				Info += ("\t" + Intensities[Pindex][Sindex]); 	
			}
			Info += ("\t" + clusterid[Pindex]);
			System.out.println(Info);
		}
	}

	public void WriteProfiles(int NumScans, String Filename, boolean OverWrite) {
		String Info;
		String Header;
		String Trailer;
		int Sindex;
		int Pindex;
		int Nindex;
		int Label;
		
                try {
                        File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);

			Header = "Scans";
                        Trailer = "RValues";
			for (Pindex = 0; Pindex < NPeaks; Pindex++) {
				Header += ("\t" + MZValues[Pindex]);
			}
			Header += "\n";
			Outgoing.write(Header);
			//write spectra before middle spectrum
			Label = -NumScans;
			for (Sindex = 1;Sindex <= NumScans; Sindex++) {
                        	Info = "" + Label;
				for (Pindex = 0; Pindex < NPeaks; Pindex++) {
                                	Info += ("\t" + Intensities[Pindex][Sindex]);
                        	}
				Info += "\n";
				Outgoing.write(Info);
				Label++;
                	}
			//write middle spectrum
                        Info = "0";
			Sindex = 0;
                        for (Pindex = 0; Pindex < NPeaks; Pindex++) {
                        	Info += ("\t" + Intensities[Pindex][Sindex]);
                        }
                        Info += "\n";
                        Outgoing.write(Info);
			//write spectra after middle spectrum
			Label = 1;
			for (Sindex = NumScans+1;Sindex <= NumScans*2; Sindex++) {
                                Info = "" + Label;
                                for (Pindex = 0; Pindex < NPeaks; Pindex++) {
                                        Info += ("\t" + Intensities[Pindex][Sindex]);
                                }
                                Info += "\n";
                                Outgoing.write(Info);
				Label++;
                        }
			Outgoing.write("\n");
			if (correlateAgainstMaster) {
				System.out.println("Writing R values");
				//write rvalues in last row
				for (Pindex = 0; Pindex < NPeaks; Pindex++) {
                        		Trailer += ("\t" + RValues[Pindex][0]);
                        	}
				Trailer += "\n";
                                Outgoing.write(Trailer);
			}
			else if (correlateAgainstAll) {
				//write rvalues as matrix
                                for (Pindex = 0; Pindex < NPeaks; Pindex++) {
					Trailer = "RValues";
					for (Nindex = 0; Nindex < NPeaks; Nindex++) {
                                        	Trailer += ("\t" + RValues[Pindex][Nindex]);
					}
					Trailer += "\n";
					Outgoing.write(Trailer);
                                }
			}
                        Outgoing.close();
                } catch (IOException failure) {
                        System.out.println("IO Error while writing " + Filename + "(WriteOutputLine)");
                        System.exit(0);
                }
        }	

        public void Correlation(int NSpectra) {
            int Pindex;
	    int Nindex;
	    //System.out.println(NPeaks);
	    if (correlateAgainstMaster) {
            	for (Pindex = 0; Pindex < NPeaks; Pindex++) {
			//System.out.println(Pindex + "\t" + NPeaks + "\t" + TopPeakIndex);
                	RValues[Pindex][0] = TwoPeakCorrelation (Pindex, TopPeakIndex, NSpectra);
			System.out.println("RValues = " + RValues[Pindex][0] + "\t" + Pindex + "\t" + TopPeakIndex);
            	}
	    }
	    else if (correlateAgainstAll) {
		for (Pindex = 0; Pindex < NPeaks; Pindex++) {
			for (Nindex = 0; Nindex < NPeaks; Nindex++) {
                        	//System.out.println(Pindex + "\t" + NPeaks + "\t" + TopPeakIndex);
                        	RValues[Pindex][Nindex] = TwoPeakCorrelation (Pindex, Nindex, NSpectra);
			}
                }
	    }
        }

        public void FindPositives() {
            int Pindex, Sindex;

            for (Pindex = 0; Pindex < NPeaks; Pindex++) {
                NPositive[Pindex] = 0;
                for (Sindex = 0; Sindex < NSpectra; Sindex++) {
                    if (Intensities[Pindex][Sindex] > 0f) {
                        NPositive[Pindex]++;
                    }
                }
            }
        }

        public void ZeroIntensity (int Sindex) {
            int Pindex;

            for (Pindex = 0; Pindex < NPeaks; ++Pindex) {
		//System.out.println(Pindex + "\t" + Sindex + "\t" + NPeaks); 
                Intensities[Pindex][Sindex] = 0f;
            }
        }
	
	public void WriteOutputLine(String Filename, String Data, boolean OverWrite) {
                try {
                        File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);
                        Outgoing.write(Data);
                        Outgoing.close();
                } catch (IOException failure) {
                        System.out.println("IO Error while writing " + Filename + "(WriteOutputLine)");
                        System.exit(0);
                }
        }

	public void WriteSpectrum (int Sindex, String Filename, boolean OverWrite, IniFile params, int Subset) {
	    int Pindex;
            int Sind = Sindex + 1;
	    String Data;
	    DecimalFormat formatter = new DecimalFormat("0.0000");
	    DecimalFormat formatter2 = new DecimalFormat("0.0");

		try {
			File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);
            
			for (Pindex = 0; Pindex < NPeaks; ++Pindex) {
				if (Subset == 1) {
                			if (RValues[Pindex][0] > params.Threshold) {
                    				Data = formatter.format(MZValues[Pindex]) + " " + formatter2.format(Intensities[Pindex][Sindex]) + /*"\t" + RValues[Pindex] + "\t" + NPositive[Pindex] +*/ "\r\n";
            					Outgoing.write(Data);
					}
				} else {
					if (RValues[Pindex][0] < params.Threshold) {
                                                Data = formatter.format(MZValues[Pindex]) + " " + formatter2.format(Intensities[Pindex][Sindex]) + /*"\t" + RValues[Pindex] + "\t" + NPositive[Pindex] +*/ "\r\n";
                                                Outgoing.write(Data);
                                        }
				}
			}
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("IO Error while writing " + Filename + "(WriteSpectrum)");
                        System.exit(0);
		}
	}	

        public void PrintSpectrum (int Sindex) {
            int Pindex;
            int Sind = Sindex + 1;

            //System.out.println("Spectrum number " + Sind);
	    System.out.println("Num Peaks Above Threshold: " + TrueNum);
            System.out.println("Num Peaks Below Threshold: " + FalseNum);
            for (Pindex = 0; Pindex < NPeaks; ++Pindex) {
//                if (RValues[Pindex] > 0.75) 
                    System.out.println(MZValues[Pindex] + "\t" + Intensities[Pindex][Sindex] + "\t\t" + RValues[Pindex][0] + "\t" + NPositive[Pindex] + "\t" + ValidatedPeak[Pindex]);
            }        
        }

	public void Percentages(IniFile params) {
	    int      	Looper;
	    
            TrueNum = 0;
	    FalseNum = 0;	
            for (Looper = 0; Looper < NPeaks; Looper++) {
                if (RValues[Looper][0] >= params.Threshold) {
                    TrueNum++;
                }
                else {
                    FalseNum++;
                }
            }
	
	}

	public void Normalizeprofiles() {
		int Sindex;
		int Pindex;
		double maxvalue = 0f;

		for (Pindex=0;Pindex<NPeaks;Pindex++) {
			maxvalue = 0f;
			//find max in signal
			for (Sindex=0; Sindex<NSpectra;Sindex++) {
				if (Intensities[Pindex][Sindex] >  maxvalue) {
					maxvalue = Intensities[Pindex][Sindex];
				}
			}
			//normalize to max
			for (Sindex=0; Sindex<NSpectra;Sindex++) {
                        	Intensities[Pindex][Sindex] = Intensities[Pindex][Sindex] / maxvalue;
			}
		}
	}

	public void Averages() {
	    int 	Looper;
	    double	TrueSum = 0f;
	    double	FalseSum = 0f;
	    int 	TrueCount = 0;
	    int 	FalseCount = 0;

	    for (Looper = 0; Looper < NPeaks; Looper++) {
		if (ValidatedPeak[Looper]) {
		    TrueCount++;
		    TrueSum += RValues[Looper][0];
		}
		else {
                    FalseCount++;
                    FalseSum += RValues[Looper][0];
                }    
	    }
	    TrueAverage = TrueSum/TrueCount;
	    FalseAverage = FalseSum/FalseCount;
	}

        public void LoadSpectrum(IniFile params, PeakList list, int sIndex) {
            	Peak	tPeak;
            	int 	pIndex;
 	    	double 	FragmentTolerance = 0.6f;
	    	boolean 	Finished = false;
	    	int 	match = 0;
		int 	cntr;
            
	    	if (sIndex == 0) {
                	System.out.println("Error: LoadSpectrum can't be called for the main spectrum!");
                	System.exit(0);
            	}
            	ZeroIntensity(sIndex); 
	    	// Daniel Cociorva 02/17/2005 - new (faster) implementation of the comparison function 
		ListIterator<Peak> peaks = list.getPeaks();
		for (pIndex = 0; pIndex < NPeaks; pIndex++) {
			Finished = false;
			match = 0;
			while(peaks.hasNext() && !Finished) {
				//System.out.println("Going Forward to Peak Index = " + peaks.nextIndex());
				tPeak = peaks.next();
				//System.out.println("First Peak in middleSpectrum = " + MZValues[pIndex]);
				//System.out.println("First Peak in currentSpectrum = " + tPeak.getM2z());
				if (tPeak.getM2z() < MZValues[pIndex] - FragmentTolerance) {
					//System.out.println("Peak MZ to Small" + "\t" + tPeak.getM2z() + "<" + MZValues[pIndex]);
				}
				else if (tPeak.getM2z() < MZValues[pIndex] + FragmentTolerance)  {
					Intensities[pIndex][sIndex] += tPeak.getIntensity();
					//System.out.println("Matched Peak" + "\t" + tPeak.getM2z() + "\t" + FragmentTolerance + "\t" + Intensities[pIndex][sIndex]);
					match++;
				}
				else {
					//System.out.println("No Match - Going Back a Few Peaks if Possible" + "\t" + tPeak.getM2z() + "\t" + MZValues[pIndex]);
					if (peaks.previousIndex() != 0) {
						for (cntr = 0; cntr < match + 1; cntr++) {
							//System.out.println("Backing Up = " + peaks.previousIndex());
							tPeak = peaks.previous();
						}
					}
					//System.out.println("Finished Going Back" + "\t" + tPeak.getM2z() + "\t" + MZValues[pIndex]);
					Finished = true;
				}
			}
		}
	}
}

class EndsWithFilter implements FilenameFilter {
        private String          FilterExtension;
        public EndsWithFilter(String extension) {
            FilterExtension = extension;
        }
        public boolean accept(File dir, String name) {
            return (name.endsWith(FilterExtension));
        }
}
