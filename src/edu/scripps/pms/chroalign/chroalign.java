package edu.scripps.pms.chroalign;

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
import edu.scripps.pms.chroalign.Elements;

public class chroalign {
	public static final int MAX_ARRAY_SIZE =1000;
	public static final int MAX_TIME_ELEMENTS = 10000;
	public static final int DEFAULTSPECTRUMSIZE = 10000;
	double[][] 		target = new double[2][MAX_TIME_ELEMENTS];
        double[][] 		reference = new double[2][MAX_TIME_ELEMENTS];
	Elements[][] 		matrix = new Elements[MAX_TIME_ELEMENTS][MAX_TIME_ELEMENTS];
	int[][]			pathArray = new int[2][2*MAX_TIME_ELEMENTS];
	int 			maxi = 0;
	int 			maxj = 0;

        public static void main(String args[]) throws Exception {
                chroalign                appObject = new chroalign();
		System.out.println("Starting Chroalign...");

		if(args.length<2) {
                        appObject.usage();
                        System.exit(0);
                }
		appObject.Go(args);
        }

	public void Go(String[] args) throws IOException {
                String          referenceMS1File;
		String          targetMS1File;
                Runtime         systemCommand = Runtime.getRuntime();
                Process         process;
	
		referenceMS1File = args[0];
		targetMS1File = args[1];
		try {
			System.out.println("Reading MS1Files ...");
			SpectrumReader refms1 = new SpectrumReader(referenceMS1File, "ms1");
			SpectrumReader tarms1 = new SpectrumReader(targetMS1File, "ms1");
	                System.out.println("Generating local distance matrix ...");
			genCostMatrix(refms1, tarms1);
			System.out.println("Generating cumulative distance matrix ...");
			genCumCostMatrix();
			System.out.println("Finding minimum path ...");
			findPath();
			System.out.println("Outputting results ...");
			writeOutputFile();

		} catch(IOException failure) {
                        System.out.println("Error while reading MS1 files");
                        System.out.println(failure.toString());
		}
        }
	
	public void genCumCostMatrix() {

		double	dmin = 0.0;
		int 	i;
		int	j;
		int	index = 0;

		for (i=1;i<maxi+1;i++) { //first for loop for rows
  			for (j=1;j<maxj+1;j++) { //second for loop for columns
				if (i!=maxi && j!=maxj) {
    					//check each of three possible points and assign minimum local cost to dmin and
					// index so can later find lowest path
					if(matrix[i][j].element <= matrix[i][j+1].element && matrix[i][j].element <= matrix[i+1][j].element) {
						index = 1;
						dmin = matrix[i][j].element;
					} else if (matrix[i][j+1].element < matrix[i][j].element && matrix[i][j+1].element < matrix[i+1][j].element) {
						index = 2;
						dmin = matrix[i][j+1].element;
					} else if (matrix[i+1][j].element < matrix[i][j].element && matrix[i+1][j].element < matrix[i][j+1].element) {
						index = 3;
						dmin = matrix[i+1][j].element;
					} else {
						System.out.println("Error calculating cumulative distances!");
						System.out.println("i " + i + " j " + j + "\t" + "1 " + matrix[i][j].element + "\t" + "2 " + matrix[i][j+1].element + "\t" + "3 " + matrix[i+1][j].element);
						System.exit(0);
					}
					matrix[i+1][j+1].element = matrix[i+1][j+1].element + dmin; //assigns cululative cost to matrix(i+1,j+1)
  					matrix[i][j].index = index;
				} else {
					matrix[i][j].index = 1;
				}
			}
		}

	
	}

	public void findPath() {
		int 	i = maxi;
		int 	j = maxj;
		int	index;
		int	cnt = 0;
		
		//start at top right hand corner of matrix and find path down to bottom left corner
		while (i > 1 && j > 1 ) {
  			index = matrix[i][j].index;
  			if (index == 1) {
    				i = i-1;
    				j = j-1;
  			} else if (index == 2) {
    				i = i-1;
  			} else if (index == 3) {
    				j = j-1;
  			} else {    
    				System.out.println("error finding path!");
				System.exit(0);
  			}
  			pathArray[0][cnt] = i;
  			pathArray[1][cnt] = j;
			cnt++;
		}	
	}

	public void genCostMatrix (SpectrumReader refms1, SpectrumReader tarms1) {
		Iterator<Peak> peaks;
		PeakList        list;
		Peak            bpeak;
		double		rt;
		int 		n;
		int 		m;
		int 		i;
		int 		j;
		boolean		useBasePeak = true;
		String[]	content = null;

		try { 
			Iterator<PeakList> refList = refms1.getSpectra();
                	Iterator<PeakList> tarList = tarms1.getSpectra();

			//put reference data into reference array
			int cnt = 0;
			while (refList.hasNext()) {
				list = refList.next();
				List<String> rtList = list.getIlines();
				content = rtList.get(0).split("\t");
				rt = new Double(content[2]).doubleValue(); 
				peaks = list.getPeaks();
				if(useBasePeak) {
					if (cnt<MAX_TIME_ELEMENTS) {
						//try to make work on base peak chromatogram first
						bpeak = basePeak(peaks);
						reference[0][cnt] = rt;
						reference[1][cnt] = bpeak.getIntensity();
					}
				}
				cnt++;
			}
			maxi = cnt;
			
			cnt = 0;
			 //put target data into target array
			while (tarList.hasNext()) {
				list = tarList.next();
				List<String> rtList = list.getIlines();
				content = rtList.get(0).split("\t");
				rt = new Double(content[2]).doubleValue();
				peaks = list.getPeaks();
				if(useBasePeak) {
					if (cnt<MAX_TIME_ELEMENTS) {
						//try to make work on base peak chromatogram first
						bpeak = basePeak(peaks);
						target[0][cnt] = rt;
						target[1][cnt] = bpeak.getIntensity();
					}
				}
				cnt++;
			}
			maxj = cnt;
		
			//create Elements objects and put into matrix
			for (i=0;i<maxi+2;i++) {
				for (j=0;j<maxj+2;j++) {
					matrix[i][j] = new Elements(); //this part is very slow should change somehow
				}
			}

			//calc euclidian distance and put into matrix
			for(i=1;i<maxi+1;i++) {
				for(j=1;j<maxj+1;j++) {
					matrix[i][j].element = euclidianDis(reference[1][i-1], target[1][j-1]);
				}
			}
			
		} catch (IOException e) {
			 System.out.println("Error while parsing MS1 files");
        	         System.out.println(e.toString());
                }
	}

	public void printMatrix() {
		int i;
		int j;
		
		for(j=0;j<matrix[1].length;j++) {
			StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
			for(i=0;i<matrix[0].length;i++) {
                                sb.append(matrix[i][j].element + "\t");
                        }
			System.out.println(sb);
                }
	}

	public void printPath() {
                int j;

                for(j=0;j<pathArray[1].length;j++) {
			if(pathArray[0][j] >= 0) {
                        	StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                        	sb.append(pathArray[0][j] + "\t" + pathArray[1][j]);
                        	System.out.println(sb);
			}
                }
        }

	public void writeOutputFile() {
                int 	j;
		int	i;
		String  TString = "0.###";
                NumberFormat    formatter = new DecimalFormat(TString);

		System.out.println("Writing output...");
		StringBuffer	xmlData = new StringBuffer();
                xmlData.append("<?xml version=\"1.0\" standalone=\"no\"?>\n");
                xmlData.append("<title>Chromatogram</title>\n");
                xmlData.append("<xLabel>Retention Time</xLabel>\n");
                xmlData.append("<yLabel>Intensity</yLabel>\n");
                xmlData.append("<noGrid/>\n");
                xmlData.append("<size width=\"550\" height=\"370\"/>\n");
                xmlData.append("<xTicks>\n");
		xmlData.append("<tick label=\"50-60\" position=\"50\" />\n");
		xmlData.append("<tick label=\"60-70\" position=\"60\" />\n");
		xmlData.append("<tick label=\"70-80\" position=\"70\" />\n");
		xmlData.append("<tick label=\"80-90\" position=\"80\" />\n");
		xmlData.append("<tick label=\"90-100\" position=\"90\" />\n");
		xmlData.append("</xTicks>\n");
		writeOutputLine("chro_out.xml",xmlData,false);
		
		StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
		System.out.println("Writing reference...");
		sb.append("<dataset name=\"reference\" marks=\"none\" connected=\"yes\" stems=\"no\">\n");
		writeOutputLine("chro_out.xml",sb,true);
		for(j=0;j<reference[1].length;j++) {
                        if(reference[1][j] > 0) {
                                sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                                sb.append("<p x=\"" + reference[0][j] + "\" y=\"" + formatter.format(reference[1][j]) + "\"/>\n");
                          	writeOutputLine("chro_out.xml",sb,true);
                        }
                }
		sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
		sb.append("</dataset>\n");
		writeOutputLine("chro_out.xml",sb,true);
		sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
		System.out.println("Writing target...");
		sb.append("<dataset name=\"target\" marks=\"none\" connected=\"yes\" stems=\"no\">\n");	
                writeOutputLine("chro_out.xml",sb,true);
		for(j=0;j<target[1].length;j++) {
                        if(target[1][j] > 0) {
                                sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                                sb.append("<p x=\"" + target[0][j] + "\" y=\"" + formatter.format(target[1][j]) + "\"/>\n");
                                writeOutputLine("chro_out.xml",sb,true);
                        }
                }
		sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                sb.append("</dataset>\n");
                writeOutputLine("chro_out.xml",sb,true);

		//Path Data
		 System.out.println("Writing path...");
		xmlData = new StringBuffer();
                xmlData.append("<?xml version=\"1.0\" standalone=\"no\"?>\n");
                xmlData.append("<title>Path</title>\n");
                xmlData.append("<xLabel>Reference Chromatogram</xLabel>\n");
                xmlData.append("<yLabel>Target Chromatogram</yLabel>\n");
                xmlData.append("<noGrid/>\n");
                xmlData.append("<size width=\"550\" height=\"550\"/>\n");
                xmlData.append("<xTicks>\n");
                xmlData.append("<tick label=\"50-60\" position=\"50\" />\n");
                xmlData.append("<tick label=\"60-70\" position=\"60\" />\n");
                xmlData.append("<tick label=\"70-80\" position=\"70\" />\n");
                xmlData.append("<tick label=\"80-90\" position=\"80\" />\n");
                xmlData.append("<tick label=\"90-100\" position=\"90\" />\n");
                xmlData.append("</xTicks>\n");
		xmlData.append("<dataset name=\"path\" marks=\"none\" connected=\"yes\" stems=\"no\">\n");
                writeOutputLine("path_out.xml",xmlData,false);
		for(j=0;j<pathArray[1].length;j++) {
                        if(pathArray[0][j] > 0) {
                                sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                                sb.append("<p x=\"" + pathArray[0][j] + "\" y=\"" + pathArray[1][j] + "\"/>\n");
                                writeOutputLine("path_out.xml",sb,true);
                        }
                }
		sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                sb.append("</dataset>\n");
                writeOutputLine("path_out.xml",sb,true);

	
        }

	public void writeOutputLine(String Filename, StringBuffer Data, boolean OverWrite) {
                try {
                        File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);
                        Outgoing.write(Data.toString());
                        Outgoing.close();
                } catch (IOException failure) {
                        System.out.println("IO Error while writing " + Filename);
                        System.exit(0);
                }
        }	

	public void printBPChro(double[][] array) {
                int i;
                int j;

                for(j=0;j<array[1].length;j++) {
			if(array[1][j] > 0) {
                        	StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
                        	sb.append(array[0][j] + "\t" + array[1][j]);
                        	System.out.println(sb);
			}
                }
        }


	public Peak basePeak(Iterator<Peak> peaks) {
		Peak	tpeak;
		Peak	bpeak = new Peak(0,0);
		double	maxint = 0.0;

		while (peaks.hasNext()) {
			tpeak = peaks.next();
			if (tpeak.getIntensity() > maxint) {
				maxint = tpeak.getIntensity();
				bpeak  = tpeak;
			}
		}
		return bpeak;
	}

	public double euclidianDis(double q, double c) {
		double d;
		d = (q - c);
		d = d*d;
		return d;
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
	
	/*private void makeMassList() {
		double tempmass = (double)startmass;
                int j = 0;
		double	increment = 0.001;
		
		while (tempmass <= endmass) {
			masslist[j] = tempmass;
			tempmass += increment;
			j++;
		}
		arraysize = j;
	}*/
        
	public void usage() {
		System.out.println("\nUSAGE: chroalign [target MS1File] [reference MS1File]");
	}
}
