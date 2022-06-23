package edu.scripps.pms.enrichmentcalc;

import java.io.*;
import java.util.*;
import java.text.*;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.SpectrumIndexReader;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.Point;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.PointList;
import edu.scripps.pms.util.spectrum.Zline;
import edu.scripps.pms.enrichmentcalc.IsotopeDistribution;

public class enrichmentcalc {
        static String             ms1Filename = null;
	public static final int	MAX_ISOTOPES = 10;
	public static final int MAX_ARRAY_SIZE = 20000;
	int	arraysize;
	int	theocnt;
	int	expcnt;
	public static final int MAX_ISOTOPE_ARRAY_SIZE = 500;
	int[]	elementcomp = new int[MAX_ISOTOPES];
	double[]         theo_abun = new double[MAX_ARRAY_SIZE];
	double[]         masslist = new double[MAX_ARRAY_SIZE];
        double[]         exp_abun = new double[MAX_ARRAY_SIZE];
 	double[] theomasses = new double[MAX_ISOTOPE_ARRAY_SIZE];
        double[] theoabun = new double[MAX_ISOTOPE_ARRAY_SIZE];
	double[] expmasses = new double[MAX_ISOTOPE_ARRAY_SIZE];
        double[] expabun = new double[MAX_ISOTOPE_ARRAY_SIZE];
	double   startmass;
        double   endmass;

        public static void main(String args[]) throws Exception {
                enrichmentcalc                appObject = new enrichmentcalc();
		int[][]		comp;
		String[]	aanames;
		
		if(args.length<4) {
                        appObject.usage();
                        System.exit(0);

                }
		appObject.Go(args);
                //appObject.Go(args, aanames, comp);
        }

	public void Go(String[] args) throws Exception {
                int             dirIterator;
		int 		i;
                String          currentMS1File;
                PeakList        list;
		Peak		tpeak;
                int             specCount = 0;
                Runtime         systemCommand = Runtime.getRuntime();
                Process         process;
		boolean		found = false;
		boolean		DTAFile = false;
		String		DTAFileName;
		boolean		OUTFile = false;
		int[]		composition;
		Iterator<Peak>	peaks;
		double	tolerance = 0.3;
		double	maxintensity = 0;
		double	enrichment;
		float	enrich_max = 100;
		float	enrich_min = 50;
		float	enrich_delta = 1;
		int	max_shift;
		int 	min_shift;
		double	threshold = 0.01;
		int	numsteps;
		double	protonMass = 1.00728;
		double	corr;
		double	resolution;
		int	scan;
		String	sequence;
		int	charge;
		double	maxcorr;
		String	mode;
		double[]	tempmass = new double[MAX_ISOTOPE_ARRAY_SIZE];
		double[]        tempint = new double[MAX_ISOTOPE_ARRAY_SIZE];

		
		numsteps = (int)((enrich_max - enrich_min) / enrich_delta);
		mode = args[0];
		currentMS1File = args[1];
		DTAFileName = args[2];
		scan = Integer.parseInt(args[3]);
		sequence = args[4];
		charge = new Integer(args[5]).intValue();
		resolution = Double.parseDouble(args[6]);
		max_shift = Integer.parseInt(args[7]);
		min_shift = -1 * Integer.parseInt(args[8]);
		enrich_max = Float.parseFloat(args[9]);
		enrich_min = Float.parseFloat(args[10]);
		enrich_delta =Float.parseFloat(args[11]);
		System.out.println("mode=" + mode + "\t" + "MS1File= " + currentMS1File + "\t" + "DTAFile= " + DTAFileName + "\t" + "sequence= " + sequence + "\t" + "scan= " + scan + "\t" + "charge= " + charge + "\t" + "resolution= " + resolution + "\t" + "maxshift= " + max_shift + "\t" + "minshift= " + min_shift + "\t" + "maxenrich= " + enrich_max + "\t" + "minenrich= " + enrich_min + "\t" + "enrichdelta= " + enrich_delta);
	
		
		//System.out.println("Reading " + currentMS1File + "...");
		SpectrumReader ms1 = new SpectrumReader(currentMS1File, "ms1");
		Iterator<PeakList> spectrumList = ms1.getSpectra();
		while (spectrumList.hasNext() && !found) {
			list = spectrumList.next();
			found = false;
			if (list.getLoscan() == scan) {
				found = true;
				System.out.println("Processing spectrum " + list.getLoscan() + "...");
				specCount++;

				//determine elemental composition
				composition = calcComposition(sequence);
				//composition = calcComposition2(sequence, aanames, comp);
				System.out.println("Elemental Composition = C" + composition[0] + " H" + composition[1] + 
					" O" + composition[2] + " N" + composition[3] + " S" + composition[4] + " P" + 
					composition[5] + " N15_" + composition[6] + " D2_" + composition[7] + 
					" C13_" + composition[8]);
			
				//call IsotopeDist to calculate start and end mass
				enrichment = enrich_max;
				IsotopeDistribution iso = new IsotopeDistribution(composition, false, enrichment/100);
				startmass = (iso.getStartMass()/charge) - 5;
				endmass = (iso.getEndMass()/charge) + 1;	
				System.out.println(startmass + "\t" + endmass);
				
				//construct mass array
				makeMassList();

				//put experimental data into arrays
				peaks = list.getPeaks();
				i = 0;
				while(peaks.hasNext()) {
					tpeak = peaks.next();
					if (tpeak.getM2z() >= startmass && tpeak.getM2z() < endmass) {
						tempmass[i] = tpeak.getM2z();
						tempint[i] = tpeak.getIntensity();
						if(tpeak.getIntensity() >= maxintensity) {
							maxintensity = tpeak.getIntensity();
						}
						i++;
					}
				}

				//normalize exp array to 1 and remove peaks below threshold
				int k = 0;
				for (int j=0;j<i;j++) {
					tempint[j] /= maxintensity;
					if (tempint[j] > threshold) {
						expabun[k] = tempint[j];
						expmasses[k] = tempmass[j];
						k++;
						expcnt = k;
					}
				}

				//transform exp data into profile data using resolution parameter
				exp_abun = generateProfile(expmasses, expabun, resolution);
				exp_abun[0] = 0;
				
				//loop through possible isotope enrichment values to find best fit
				System.out.println("Starting Correlation Analysis");
				System.out.println("% Enrichment\tCorrelation\tShift Required");	
				while(enrichment >= enrich_min) {
					iso = new IsotopeDistribution(composition, false, enrichment/100);
					tempmass = iso.getMasslist();
					tempint = iso.getRelAbun();
					int j = 0;
					for(i=0;i<tempmass.length;i++) {
						if (tempint[i] > threshold) {
							theomasses[j] = (tempmass[i] + (charge * protonMass)) / 
								charge;
							theoabun[j] = tempint[i];
							theocnt = j;
							j++;
						}
					}
					
					//transform into profile data using resolution parameter
					theo_abun = generateProfile(theomasses, theoabun, resolution);
					
					//for (i=0;i<arraysize;i++) {
                        		//	System.out.println(masslist[i] + "\t" + exp_abun[i] + "\t" + theo_abun[i]);
                			//}
					//correlate theoretical and experimental data
					maxcorr = -100;
					int shiftmax = 0;
					for (int p = min_shift; p<max_shift+1;p++) {
						corr = correlation(theo_abun,exp_abun,arraysize, p);
						//corr = correlation2(theomasses,theoabun,expmasses,expabun,MAX_ISOTOPE_ARRAY_SIZE, tolerance);
						if(corr>maxcorr) {
							maxcorr = corr;
							shiftmax = p;
						}
					}
					
					//output correlation results
					System.out.println(enrichment + "\t" + maxcorr + "\t" + shiftmax);
					
					//decrement enrichment
                                        enrichment -= enrich_delta;
				}	
			}

		}
        }

	public double correlation (double[] abun1, double[] abun2, int number, int shift) {
            	double r = 0;
            	int index;
            	double Sumx, Sumy, Sumxy, Sumx2, Sumy2, x, y;

            	Sumx = 0; Sumy = 0; Sumx2 = 0; Sumy2 = 0; Sumxy = 0;
            	for (index = 0; index < number; index++) {
			if (index + shift > 0 && index + shift < number) {
                		x = abun1[index + shift];
			} else {
				x = abun1[index];
			}
                	y = abun2[index];
                	Sumx += x;
                	Sumy += y;
                	Sumx2 += x*x;
                	Sumy2 += y*y;
                	Sumxy += x*y;
            	}
            	r = (Sumxy - Sumx*Sumy/number)/((double) Math.sqrt((Sumx2 - Sumx*Sumx/number) * (Sumy2 - Sumy*Sumy/number)));
		return r;
        }
	
	private double[] generateProfile(double[] mass, double[] intensity, double resolution) {
		double		deltam;
		double[]	result = new double[MAX_ARRAY_SIZE];
		int 		j = 0;
		int 		i = 0;
		double		secondterm;
		double		maxresult = 0;
		int		count = 0;
		while (mass[j] > 0) {
			count++;
			j++;
		}
		
		//determine peak width from resolution
		for (j=0;j<count;j++){
			
			deltam= mass[j]/resolution;
			i=0;
			while (i<arraysize) {
				secondterm = 0;
				secondterm = ((masslist[i] - mass[j])*(masslist[i] - mass[j])) / (2 * deltam*deltam);
				secondterm = Math.pow(2.71828,(-1 * secondterm));
				secondterm /= deltam * Math.sqrt((2*3.14159));
				secondterm = secondterm * intensity[j];
				result[i] = secondterm + result[i];
				if (result[i] > maxresult) {
					maxresult = result[i];
				} 
				i++;
		       }
            	}
		// normalize data to 1
		for (i=0;i<arraysize;i++) {
			result[i] /= maxresult;
			//System.out.println(masslist[i] + "\t" + result[i]);
		}
		return result;
	}

	private void makeMassList() {
		double tempmass = (double)startmass;
                int j = 0;
		double	increment = 0.001;
		
		while (tempmass <= endmass) {
			masslist[j] = tempmass;
			tempmass += increment;
			j++;
		}
		arraysize = j;
	}
	
	
	public int[]  calcComposition2 (String sequence, String[] aanames, int[][] comp) {

		//Calculate the number of each amino acid present in the sequence  
		int	a;
		int	i;
		aaComp	sequenceComp = new aaComp();
		List<AA>	aaList = new ArrayList();
		AA	aa;
		
		//aaList = sequenceComp.readAAComp();
		//for each character in sequence aum elemental comp
		for(i=0;i<sequence.length();i++) {
			for (int j=0;j<aanames.length;j++) {
				if(sequence.substring(i,i+1).equals(aanames[j])) {
					sumElements2(comp, j);
				}
			}
		}
		//add 2 protons and one oxygen for termini
		elementcomp[1]+=2;
		elementcomp[2]++;
			
		return elementcomp;
	}

	public int[]  calcComposition (String sequence) {
                String comp;

                //Calculate the number of each amino acid present in the sequence
                int     a;
                int     i;
                aaComp  sequenceComp = new aaComp();
                List<AA>        aaList = new ArrayList();
                AA      aa;

                aaList = sequenceComp.readAAComp();
                //for each character in sequence aum elemental comp
                for(i=0;i<sequence.length();i++) {
                        for (int j=0;j<aaList.size();j++) {
                                aa = aaList.get(j);
                                if(sequence.substring(i,i+1).equals(aa.aa)) {
                                        sumElements(aa.elements);
                                }
                        }
                }
                //add 2 protons and one oxygen for termini
                elementcomp[1]+=2;
                elementcomp[2]++;

                return elementcomp;
        }

        private void sumElements(int[] arr){
                for(int i=0; i<MAX_ISOTOPES; i++)
                {
                    elementcomp[i] += arr[i];
                }
        }		

	private void sumElements2(int[][] arr, int j){
		for(int i=0; i<MAX_ISOTOPES; i++)
		{
		    elementcomp[i] += arr[j][i];
		}
    	}
 
	public void usage() {
		System.out.println("\nUSAGE: enrichmentcalc [MS1File] [Scan] [Sequence] [Charge] [flags(optional)]");
  		System.out.println("\nFlags must be individually called with a dash. Example: -f");
  		System.out.println("Valid flags are:");
  		System.out.println("\t-f [filename], reads DTAselect-filter.txt file as input\n");
  		System.out.println("\t-o [filename], outputs results to file\n");
	}
}

class AA {
	private static final int	MAX_ISOTOPES = 10;
	String	aa;
	int[]	elements = new int[MAX_ISOTOPES];
}

class aaComp {
	int	i;
	List<AA>	aaList = new ArrayList();
	private static final int        MAX_ISOTOPES = 10;
	
        public List<AA> readAAComp() {
                try {
                        //defines local variables and checks to see that the precursorcalc.params file exists
                        File            comp = new File("AAComposition.ini");
			//System.out.println("Reading AAComposition.ini");
                        if (comp.exists()) {
                                FileReader      InputFileReader = new FileReader(comp);
                                BufferedReader  Incoming = new BufferedReader(InputFileReader);
                                String          LineBuffer;
				String		tempString;
                                String          WholeLine;
                                StringTokenizer Parser;
                                WholeLine = Incoming.readLine();
                                //Reads in each line and assigns values to appropriate instance variables
                                while (WholeLine != null) {
                                        Parser = new StringTokenizer(WholeLine, "\t");
                                        if (Parser.hasMoreTokens()) {
                                                LineBuffer = Parser.nextToken();
                                                if (LineBuffer.startsWith("#") || LineBuffer.startsWith("AA")) {
                                                        // It's a comment; ignore it.
                                                }
                                                else if (LineBuffer.matches("<[A-Z]>")) {
							AA	tempAA = new AA();
							Parser = new StringTokenizer(WholeLine, "\t");
							tempString = Parser.nextToken();
							for (i=0;i<MAX_ISOTOPES;i++){
								addAA(tempAA, LineBuffer.substring(1,2), new Integer(Parser.nextToken()).intValue());
							}
							aaList.add(tempAA);
                                                }	
                                                else {
                                                        System.out.println("I don't understand this option in AAComposition.ini.");
                                                        System.out.println(WholeLine);
                                                        System.exit(0);
                                                }
                                        }
                                        WholeLine = Incoming.readLine();
                                }
                                // Close file
                                Incoming.close();
                        } else {
				System.out.println("AAComposition.ini file not found");
			}
                }
                //error handling
                catch (IOException failure) {
                        System.out.println("Error while reading AAComposition.ini");
                        System.out.println(failure.toString());
                }
                //return this;
		return aaList;
        }

	public void addAA(AA tempAA, String name, int value) {
		
		tempAA.aa = name;
		tempAA.elements[i] = value;
	}
	
        public void printComp(List<AA> list) {
		AA	aa;
                System.out.println("Amino Acid Composition");
		System.out.println("Composition File has " + list.size() + " entries");
		for (int i = 0; i<list.size();i++) {
			aa = list.get(i);
			System.out.println("<" + aa.aa + ">");
			for (int j = 0;j<MAX_ISOTOPES;j++) {
				System.out.println(aa.elements[j]);
			}
		}	
        }
	
	
}
