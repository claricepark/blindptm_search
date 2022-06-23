package edu.scripps.pms.qcorr;

import java.io.*;
import java.util.*;
import java.text.*;

/*
    program for testing the cross correlation between a spectrum and a reversed copy of the same spectrum.
    John Venable, december 9, 2004.
*/
public class QCorr {
    	static String             ms2Filename = null;
	String          	  homedir = System.getenv("QCORR_HOME");
    	public ParamsFile         params = new ParamsFile(homedir);
    	
	public static void main(String args[]) throws Exception {
        	QCorr appObject = new QCorr();
		appObject.Go(args);
    	}

	public void Go(String[] args) throws Exception {
		File            currentDirectory = new File(System.getProperty("user.dir"));
                EndsWithFilter  uniMS2Filter = new EndsWithFilter(".ms2");
		EndsWithFilter  uniMS1Filter = new EndsWithFilter(".ms1");
                String[]        uniMS2List = currentDirectory.list(uniMS2Filter);
		String[]        uniMS1List = currentDirectory.list(uniMS1Filter);
                int             dirIterator;
                String          currentMS2File;
		String          currentMS1File;
		PeakList        list;
		String[]	data;
		ArrayList<String>	spectrumData;
		int		specCount = 1;
		boolean		filter = false;
		int		minPeaks = 0;
		Runtime		systemCommand = Runtime.getRuntime();
		Process		process;
		SimplifySpectra	simplify = new SimplifySpectra();
		boolean		profileMode;
		List<String>	hLines = new ArrayList();
		SpectrumReader ms1;
		Iterator<PeakList>	ms1SpectraList;


		 //ReadParams method
                //params.readParams(homedir);
		params.printParams();
                
		//read in MS2 and Parse through MS2. Each spectrum is processed in turn
		for(dirIterator = 0; dirIterator < args.length; dirIterator++) {
			currentMS2File = args[dirIterator];
			System.out.println("Reading " + currentMS2File + "...");
			SpectrumReader ms2 = new SpectrumReader(currentMS2File, "ms2");
			Iterator<PeakList> spectrumList = ms2.getSpectra();
			ms1SpectraList = ms2.getSpectra();
			if (params.AssignMWFromMS1) {
				if (dirIterator < uniMS1List.length) {
					currentMS1File = uniMS1List[dirIterator];
					System.out.println("Reading " + currentMS1File + "...");
					ms1 = new SpectrumReader(currentMS1File, "ms1");
					ms1SpectraList = ms1.getSpectra();
				} else {
					System.out.println("I can't find appropriate MS1 file");
					System.exit(0);
				}
			}
			Spectrum2 spectrum = new Spectrum2();
			while (spectrumList.hasNext()) {
				list = spectrumList.next();
				hLines = list.getHlines();
				System.out.print("Processing spectrum " + list.getLoscan() + "...\r");
				//if(list.getLoscan()==15420) {
				//	System.out.println(list.getLoscan() + "\t" + filter);
				//}
				if (list.numPeaks() > minPeaks) {
					if (params.DeChargeSpectrum) {
						//is data profile or centroided
						profileMode = dataMode(list);
						list = simplify.simplifySpectra(list, params.IsotopeAccuracy, profileMode);
					}
					filter = false;
					data = spectrum.processSpectrum2(params, list, specCount, currentMS2File, ms1SpectraList);
					//System.out.println("\nSpectra processed so far:\t" + specCount);
					filter = createSpectrumData(params, list, specCount, data);
					if (!filter) {
						spectrum.writeSpectrum(list, currentMS2File, specCount, hLines);
						if (params.OutputFile) {
							 spectrum.writeData(data, currentMS2File + "_Precursors.txt", list, specCount);
						}
					}
					//if (params.OutputFile) {
					//	spectrum.writeData(data, currentMS2File + "_Precursors.txt", list, specCount);
					//}
				}
				specCount++;

			}
			if (params.OverWriteMS2) {
				String tempString = "mv " + "QCorr_" + currentMS2File + " " + currentMS2File;
				System.out.println(tempString);
				process = systemCommand.exec("mv " + "QCorr_" + currentMS2File + " " + currentMS2File);
			}
		}
	}

	boolean dataMode(PeakList list) {
		boolean mode = false;
		int 	i;
		Peak	p = new Peak(0f,0f);
		float	mz1;
		float	mz2;
		ListIterator<Peak> peaks = list.getPeaks();
		
		p = peaks.next();
		mz1 = p.getM2z();
		p = peaks.next();
                mz2 = p.getM2z();
		if (mz2 - mz1 < 0.1) {
			mode = true;
		}
		return mode;		
	}

	boolean createSpectrumData(ParamsFile params, PeakList list, int specCount, String[] data) {
		ArrayList<String> 	spectrumList = new ArrayList();
		String			tempSpectrum;
		StringTokenizer		parser;
		String			lineBuffer;
		int			i;
		int			cnt;
		float[]			molWeight = new float[data.length];
		float[]                 qCorr = new float[data.length];
		float[]                 mzs = new float[data.length];
		int[]                 charges = new int[data.length];
		ArrayList<Zline>	newZlines = new ArrayList();
		Zline			z;
		Iterator<Zline>		tZLine;
		String			dLine;
		boolean			filter = false; 
		float			threshPlus1 = 0;
		float                   threshPlus2 = 0;
		float                   threshPlus3 = 0;
		int			chg = 0;
		float			tempfloat = 0;
		int			thresh = 0;
		
		for (cnt = 0; cnt < data.length; cnt++) {
			parser = new StringTokenizer(data[cnt],"\t");
			molWeight[cnt] = new Float(parser.nextToken()).floatValue();
			qCorr[cnt] = new Float(parser.nextToken()).floatValue();
		}	
		Iterator<Zline> zLines = list.getZlines();
		if(params.FilterSpectra) {
			newZlines = new ArrayList();	
			//remove zLine if QCorr < params.Threshold
			for (cnt = 0; cnt < data.length; cnt++) {
				z = zLines.next();
				chg = z.getChargeState();
				if(chg==1) {
					thresh = 0;
				} else if (chg==2) {
					thresh = 1;
				} else if (chg==3) {
					thresh = 2;
				} else {
					thresh = 3;
				}
				if(qCorr[cnt] > params.FilterThreshold[thresh]) {;
					if (params.ReplaceMW) {
						tempfloat = molWeight[cnt] + 1.00728f;
                                        	z.setM2z(tempfloat);
                                	}
					dLine = "D\tQCorr\t" + qCorr[cnt];
					z.addDline(dLine);
                                	dLine = "D\tMolWeight\t" + molWeight[cnt];
                                	z.addDline(dLine);
					newZlines.add(z);		
				}
			}
			list.setZlines(newZlines);
			/*tZLine = list.getZlines();
			while(tZLine.hasNext()){
				z = tZLine.next();
				System.out.println(z.getChargeState() + "\t" + z.getM2z());
			}*/
		} else {
			//add QCorr score and molecular weight to D lines
			for (cnt = 0; cnt < data.length; cnt++) {
				z = zLines.next();
				//replace molecular weight
				if (params.ReplaceMW) {
					tempfloat = molWeight[cnt] + 1.00728f;
					z.setM2z(tempfloat);
				}
				dLine = "D\tQCorr\t" + qCorr[cnt];
                        	z.addDline(dLine);
				dLine = "D\tMolWeight\t" + molWeight[cnt];
                                z.addDline(dLine);
				newZlines.add(z);
                        }
			list.setZlines(newZlines);
		}
		//if all charge states score below thresholds then remove whole spectrum
		if (newZlines.size() == 0) {
                	filter = true;
                } else {
			if(params.AssignCharge) {
				//keeps only highest scoring charge state by removing lowest scoring Z lines
				//find highest scoring charge state 
				float max = qCorr[0];
				int maxcnt = 0;
				newZlines = new ArrayList();
				for (cnt = 0; cnt < data.length; cnt++) {
                                	if(qCorr[cnt] > max) {
						max = qCorr[cnt];
                                        	maxcnt = cnt;
                                	}
                        	}
				zLines = list.getZlines();
				for (cnt = 0; cnt < data.length; cnt++) {
					z = zLines.next();
					if (cnt == maxcnt) {
						newZlines.add(z);
                        			list.setZlines(newZlines);
					}	
				}
			}
		}
		
		return filter;
	}
}

class Spectrum2 {
	String[] processSpectrum2(ParamsFile params, PeakList list, int specCount, String currentMS2File, Iterator<PeakList> ms1SpectraList) {
		ArrayList<Peak>	originalSpectrum;
		ArrayList<Peak> reversedSpectrum;
		float[]	correlationResults;
		int jj = 0;
		String[] tempData = new String[1];
                StringBuffer directory = new StringBuffer();
		int	chgcnt = 0;
		Zline tempZline;
		File f = new File(".");
		int	ARRAY_NUM;
		int	p = 0;
		int	num;
		ArrayList<Point> molWeights = new ArrayList();
		PointList	cenResults;
		double refinedMZ;
		String 	scan;
		SpectrumIndexReader reader;
		String  fileName;
		PointList[] ccResults;
		Point pTemp;
		ListIterator<Point> pList;
                int     tao;
		List<Point> cenScores;

		
                fileName = currentMS2File.substring(0,currentMS2File.length() - 3) + "ms1";
		//Determine chargeStates
                Iterator<Zline> zLines = list.getZlines();
                ArrayList<Integer> charges = new ArrayList<Integer>();
                while(zLines.hasNext()) {
                        tempZline = zLines.next();
                        charges.add(new Integer(tempZline.getChargeState()) );
                        chgcnt++;
                }
		directory.append(f.getAbsolutePath());
		
		//calculate ccResults using cross correlation
		ccResults = list.calcQCorrs(params.PeakThreshold, (int)(1/params.MassAccuracy) , (double)params.InjWindow);

		//calculate MW for each charge state
		tempData = new String[chgcnt];
		for (jj = 0; jj < chgcnt; jj++) {
			//add smallest potential MW and largest potential MW to mzData
			ArrayList mzData = new ArrayList();
			pList = ccResults[jj].getPoints();
			pTemp = pList.next();
			mzData.add(pTemp.getXValue());
			while(pList.hasNext()) {
                        	pTemp = pList.next();
			}
                        mzData.add(pTemp.getXValue());
			//centroid ccResults
			cenResults = centroidCCArray(ccResults[jj], mzData, params);
			//write CCearray to text output file if OutputArrays requested
			if (params.OutputCCArray) {
				printCCScores(list, mzData, ccResults[jj], specCount, currentMS2File, directory, 
					charges.get(jj).intValue(), params, cenResults);
			}
			cenScores = cenResults.getSortedPoints(true);
			molWeights = finalizePrecursorSelection(params, cenScores);
			//if using high mass accuracy full scan to assign final MW then use this routine
			if (params.AssignMWFromMS1) {
				refinedMZ = (molWeights.get(0).getXValue() + (charges.get(jj).intValue() * 1.00728)) 
					/charges.get(jj).intValue();
				//determine appropriate MS spectrum number (i.e., MS spectrum immediately preceeding identification)
				scan = determineMSScan(fileName, list.getLoscan());
				//get MS1 spectra info using .index file
				try {
					reader = new SpectrumIndexReader(fileName, scan);
					//generate PointList containing MS data
					PointList msResults = reader.getPointList(); 
					//call function (like centroidCCArray) to assign and centroid peaks
					PointList msPeaks;
					msPeaks = centroidMSData(msResults, molWeights.get(0).getXValue(), params, 0f, 
						charges.get(jj).intValue());
					//assign new MW values to molWeights
                        		cenScores = msPeaks.getSortedPoints(true);
					tempData[jj] = cenScores.get(cenScores.size()-1).getXValue() + "\t" + 
						cenScores.get(cenScores.size()-1).getIntensity();
				} catch (IOException e) {
                        		System.out.println("index file does not exist " + e);
                		        System.exit(0);
		                }
			} else {
				tempData[jj] = molWeights.get(0).getXValue() + "\t" + molWeights.get(0).getIntensity();
			}
		}
		return tempData;
	}

	String determineMSScan(String fileName, int scan) {
		PeakList	ms1List;
		String 		max = "";
		PeakList	list;
		boolean	found = false;
		String lastLine;
		String[]	strArr = new String[2];
        	StringBuffer sb = new StringBuffer();
		BufferedReader br;
		
		try {
        		br = new BufferedReader(new FileReader(fileName + ".index"), 4096);
                	while ((lastLine = br.readLine()) != null && !found) {
                        	strArr = lastLine.split("\t");
                        	if (Integer.parseInt(strArr[0]) <= scan) {
                                	max = strArr[0];
                        	} else {
                                	found = true;
                        	}
                	}
                	if(null == lastLine) {
                        	System.out.println("Spectrum number does not exist : ");
                	}
		} catch(IOException e) {
                	System.out.println("index file does not exist " + e);
                	//throw new IOException("index file does not exist " + e);
			System.exit(0);
        	} 
		return max;
	}

	PointList centroidMSData(PointList msResults, double MW, ParamsFile params, double peakThreshold, int charge) {
                double cenThreshold = 0.3f;
                double offset = cenThreshold;
                double avgIntensity;
                double sumxy;
                double sumy;
                int     i;
                double xbar;
                List<PointList> ccPeaks;
                PointList       pList;
                Point   p;
		double accuracy = 0.3f;
                Iterator<Point> points;
                PointList tempList = new PointList();
                PointList tempList2;
		tempList2 = findRoughMZ(msResults, MW ,params, charge, accuracy);
		ccPeaks = findPeaks(tempList2, peakThreshold);
		Point p2 = new Point(MW - accuracy ,0);
                tempList.addPoint(p2);
                for(i=0;i<ccPeaks.size();i++) {
                        pList = ccPeaks.get(i);
                        points = pList.getPoints();
                        sumxy = 0;
                        sumy = 0;
			System.out.println("new peak");
                        while(points.hasNext()) {
                                p = points.next();
				System.out.println("centroid points" + p.getXValue() + p.getIntensity());
                                if ( p.getIntensity() > (cenThreshold * pList.getMaxIntensity()) ) {
                                        sumxy += (p.getXValue() * p.getIntensity());
                                        sumy += p.getIntensity();
					System.out.println("Point added to centroid"); 
                                }
                        }
                        xbar = sumxy/sumy;
			System.out.println("centroid = " + xbar);
                        p2 = new Point(xbar,sumy);
                        tempList.addPoint(p2);
                }
                p2 = new Point(MW + accuracy,0);
                tempList.addPoint(p2);
		return tempList;
	}		

	PointList centroidCCArray(PointList ccResults, ArrayList mzData, ParamsFile params) {
		double peakThreshold;
		double cenThreshold = 0.3f;
		double offset = cenThreshold;
		double avgIntensity;
		double sumxy;
		double sumy;
		int	i;
		double xbar;
		List<PointList>	ccPeaks;
		PointList	pList;
		Point   p;
		Iterator<Point> points;
		PointList tempList = new PointList();
		PointList tempList2 = new PointList();
		//System.out.println("Starting Centroiding Routine");
		avgIntensity = getAverage(ccResults);
		//System.out.println("Avg Intensity = " + avgIntensity);
		peakThreshold = offset * ccResults.getMaxIntensity();
		//peakThreshold = offset * avgIntensity;
		//System.out.println("Threshold = " + peakThreshold);
		ccPeaks = findPeaks(ccResults, peakThreshold);
		//System.out.println("Number of Peaks Found = " + ccPeaks.size());
		Point p2 = new Point(((Double)mzData.get(1)).doubleValue(),100);
                tempList2.addPoint(p2);
		for(i=0;i<ccPeaks.size();i++) {
			pList = ccPeaks.get(i);
			points = pList.getPoints();
			sumxy = 0;
			sumy = 0;
			while(points.hasNext()) {
				p = points.next();
				if ( p.getIntensity() > (cenThreshold * pList.getMaxIntensity()) ) {
					sumxy += (p.getXValue() * p.getIntensity());
					sumy += p.getIntensity();
				}
			}
			xbar = sumxy/sumy;
			p2 = new Point(xbar,sumy);
			tempList2.addPoint(p2);
		}
		p2 = new Point(((Double)mzData.get(0)).doubleValue(),100);
                tempList2.addPoint(p2);
		return tempList2;
	}

	PointList findRoughMZ(PointList msResults, double MW, ParamsFile params, int charge, double accuracy) {
		Point   p;
                Point   p2;
                Iterator<Point> msData = msResults.getPoints();
                PointList       tempList = new PointList();
		double	pMW;

		while(msData.hasNext()) {
                        p = msData.next();
			pMW = (p.getXValue() * charge - (charge * 1.00728));
			if (pMW > MW - accuracy && pMW <= MW + accuracy) {
				System.out.println("Peak Range" + p.getXValue() + "\t" + pMW);
				p.setXValue(((p.getXValue()*charge)-(charge*1.00728)));
				tempList.addPoint(p);
			}
                }
                return tempList;
	}

	List<PointList> findPeaks(PointList ccResults, double threshold) {
		Point   p;
		List<PointList>	list = new ArrayList();
		boolean	newPeak = true;
		PointList points = new PointList();
		
		Iterator<Point> ccData = ccResults.getPoints();
		while(ccData.hasNext()) {
                        p = ccData.next();
                        if (p.getIntensity() > threshold) {
				if(newPeak) {
					points = new PointList();
				}
				newPeak = false;
				while (p.getIntensity() > threshold && ccData.hasNext()) {
					points.addPoint(p);
					p = ccData.next();
				}
				list.add(points);
			} else {
				newPeak = true;
			}
                }
		return list;
	}

	float getAverage(PointList ccResults) {
		float   avgIntensity;
		float	sum = 0;
		Point	p;

		Iterator<Point> ccData = ccResults.getPoints();
		while(ccData.hasNext()) {
			p = ccData.next();
			sum += p.getIntensity();
		}
		avgIntensity = sum / ccResults.numPoints();
		return avgIntensity;
	}

	public void writeSpectrum(PeakList list, String currentMS2File, int specCount, List<String> hLines) {
		String	tempSpectrum;
		String	tempHlines = new String();
		int	i;
			
                tempSpectrum = list.getSpectrumWithoutHlines();
		for (i=0;i < hLines.size();i++) {
			tempHlines = tempHlines + hLines.get(i) + "\n";
		}
		if(specCount == 1) {
			tempSpectrum = tempHlines + tempSpectrum;
			writeDataLine("QCorr_" + currentMS2File, tempSpectrum, false);
                } else {
                        writeDataLine("QCorr_" + currentMS2File, tempSpectrum, true);
                }
	}
	
	ArrayList finalizePrecursorSelection(ParamsFile params, List<Point> ccScores) {
		double answer;
		double tempanswer;
		int	arrcntr;
		double finalMZ = 0f;
		double finalIntensity = 0f;
		ArrayList<Point> 	tempList = new ArrayList();
		Point	p;
			
		//remove 2 protons (+1 b ion X +1 Y ion equals M + 2) 
		finalMZ = ccScores.get(ccScores.size() - 1).getXValue() - 2.01456;
		answer = finalMZ;
		finalIntensity = ccScores.get(ccScores.size() - 1).getIntensity();
		
		//check remaining peaks to see if monoisotopic peak is present
		for (arrcntr = 1; arrcntr < params.NumIons; arrcntr++) {
			if ((ccScores.size() - 1 - arrcntr) > 0) {
				tempanswer = ccScores.get(ccScores.size() - 1 - arrcntr).getXValue() - (float)2.01456;
				//see if peak is present within mass accuracy tolerance
				if ((answer - tempanswer) < 1.00728 + 0.1 && (answer - tempanswer) > 1.00728 - 0.1) {
					//see if peak is > 50% of max intensity
					if (ccScores.get(ccScores.size() - 1 - arrcntr).getIntensity() >= (ccScores.get(ccScores.size() - 1).getIntensity() * 0.5)) {
						//System.out.println("isotope peak found!\n");
						finalMZ = tempanswer;
						finalIntensity = ccScores.get(ccScores.size() - 1 - arrcntr).getIntensity();
					} 
				}
			}
		}
		//System.out.println("monoisotopic mass = " + finalanswer + "\n");
		p = new Point(finalMZ, finalIntensity);
		tempList.add(p);
		return tempList;	
	}

	public void printSpectrum(ArrayList<Peak> spectrum) {
                int     cnt;
		float	tempmz;
		float	tempint;

                for (cnt = 0; cnt < spectrum.size(); cnt++) {
			tempmz = spectrum.get(cnt).getM2z();
			tempint = spectrum.get(cnt).getIntensity();
                	System.out.println(tempmz + "\t" + tempint);
		}
	}

	public void printCCSpectrum (ListIterator<Point> spectrum, String fileName, ArrayList mzData, ParamsFile params, ListIterator<Point> spectrum2) {
		Point	tempPoint;
		Point	tempPoint2;
		String	data;
		String	temp;
		String  tempint;

		double	finalMZ = ((Double)mzData.get(0)).doubleValue();

		while (spectrum.hasNext()) {
                        tempPoint = spectrum.next();
			temp = FormatFloatNum(new Float(finalMZ),5);
			tempint = FormatFloatNum(new Float(tempPoint.getIntensity()), 1);
                        data = temp + "\t" + tempint + /*"\t" + tempPoint.getXValue() */"\r\n";
			writeDataLine(fileName + "_CCArrays.txt", data, true);
			finalMZ -= params.MassAccuracy;
                }
		while (spectrum2.hasNext()) {
                        tempPoint2 = spectrum2.next();
                        //data = tempPoint2.getXValue() + "\t" + tempPoint2.getIntensity() + "\r\n";
                        //writeDataLine(fileName + "_centroided.txt", data, true);
                }
	}

	public void printCCScores (PeakList list, ArrayList mzData, PointList ccResults, int specCount, String fileName, 
		StringBuffer directory, int chargeState, ParamsFile params, PointList cenResults) {
               	
		Point ptemp;

                ListIterator<Point> pList = ccResults.getPoints();
		ListIterator<Point> p2List = cenResults.getPoints();
		if (specCount == 1 && chargeState < 3) {
			//write CC array to text file
		    	writeDataLine(fileName + "_CCArrays.txt", "S" + "\t" + list.getLoscan() + "\t"+ chargeState + "\t" + "\r\n", false);

			//write Centroided array to text file
			/* writeDataLine(fileName + "_centroided.txt", "S" + "\t" + list.getLoscan() + "_" + chargeState + "\t" 
				+ list.getLoscan() + "\t"+ chargeState + "\t" + "\r\n", false);

	    		//write html output file if OutputArrays requested
			writeDataLine(fileName + "_CCArrays.html", "<a href=\"http://bart.scripps.edu/pages/SpectrumViewer/?fileName=" + 
				directory + "/" + fileName + "_CCArrays.txt" + "&spectrumNum=" + list.getLoscan() + "_" + 
				chargeState + "\">" + list.getLoscan() + "_" + chargeState + "</a><br>\r\n",false);

			//write html output file if OutputArrays requested
			writeDataLine(fileName + "_centroided.html", "<a href=\"http://bart.scripps.edu/pages/SpectrumViewer/?fileName=" 
				+ directory + "/" + fileName + "_centroided.txt" + "&spectrumNum=" + list.getLoscan() + "_" +
                        	chargeState + "\">" + list.getLoscan() + "_" + chargeState + "</a><br>\r\n",false);*/
		} else {
			//write CC array to text file
			writeDataLine(fileName + "_CCArrays.txt", "S" + "\t" + list.getLoscan() + "\t"+ chargeState + "\t" + "\r\n", true);

			//write Centroided array to text file
                        /* writeDataLine(fileName + "_centroided.txt", "S" + "\t" + list.getLoscan() + "_" + chargeState + "\t" + 
				list.getLoscan() + "\t"+ chargeState + "\t" + "\r\n", true);

			//write html output file if OutputArrays requested
			writeDataLine(fileName + "_CCArrays.html", "<a href=\"http://bart.scripps.edu/pages/SpectrumViewer/?fileName=" +
                        	directory + "/" + fileName + "_CCArrays.txt" + "&spectrumNum=" + list.getLoscan() + "_" +
                        	chargeState + "\">" + list.getLoscan() + "_" + chargeState + "</a><br>\r\n",true);

			//write html output file if OutputArrays requested
                        writeDataLine(fileName + "_centroided.html", "<a href=\"http://bart.scripps.edu/pages/SpectrumViewer/?fileName=" + 
				directory + "/" + fileName + "_centroided.txt" + "&spectrumNum=" + list.getLoscan() + "_" +
                        	chargeState + "\">" + list.getLoscan() + "_" + chargeState + "</a><br>\r\n",true);*/

		} 
		//Print each point in ccList
		printCCSpectrum(pList, fileName, mzData, params, p2List); 
        }
	
	public void writeData(String[] data, String fileName, PeakList list, int specCount) {
                int dataloop;
		Iterator<Zline> zLines;
		Zline	z;
		String	calcMW;
		String	CCscore;
		String[]        strArr = new String[2];	

		zLines = list.getZlines();
                for (dataloop = 0; dataloop < data.length; dataloop++) {
			if (zLines.hasNext()) {
                        	//newZlines = new ArrayList();
                        	z = zLines.next();
                        	if (specCount == 1 && dataloop == 0) {
					writeOutputLine(fileName, "Scan\tCharge\tCalculated MW\tCC Score" + "\n", false);
					strArr = data[0].split("\t");
					calcMW = FormatFloatNum(new Float(strArr[0]),5);
					CCscore = FormatFloatNum(new Float(strArr[1]),1);
                                	data[0] = list.getLoscan() + "\t" + z.getChargeState() + "\t" + calcMW + "\t" + CCscore;
                                	writeOutputLine(fileName, data[0] + "\n", true);
                        	}
                        	else {
					strArr = data[dataloop].split("\t");
					calcMW = FormatFloatNum(new Float(strArr[0]),5);
                                	CCscore = FormatFloatNum(new Float(strArr[1]),1);
                                	data[dataloop] = list.getLoscan() + "\t" + z.getChargeState() + "\t" + calcMW + "\t" + CCscore;
                                	writeOutputLine(fileName, data[dataloop] + "\n", true);
                        	}
                	}
		}
        }

        public void writeDataLine(String fileName, String data, boolean keepOrig) {

                if (!keepOrig) {
                        writeOutputLine(fileName, data, false);
                }
                else {
                        writeOutputLine(fileName, data, true);
                }
        }

        public void writeOutputLine(String Filename, String Data, boolean OverWrite) {
                try {
                        File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);
                        Outgoing.write(Data);
                        Outgoing.close();
                } catch (IOException failure) {
                        System.out.println("IO Error while writing " + Filename);
                        System.exit(0);
                }
        }	

	public void Print1DArray (float[] Array) {
                int cnt = 0;

                while (cnt < Array.length) {
                        System.out.println(Array[cnt]);
                        cnt++;
                }
        }

	String FormatFloatNum(float number, int decimalplaces) {
		String	TString = "0.";
		int	i;

		for (i = 0;i < decimalplaces;i++) {
			TString = TString + "#";
		}
		NumberFormat	formatter = new DecimalFormat(TString);
		TString = formatter.format(number);
		return TString;
	}
	
	ArrayList<Peak> normalizeSpectrum(ParamsFile params, ArrayList<Peak> binnedSpectrum, PeakList list) {
		float	max = 0;
		int looper = 0;
		float	normInt;
		float	mz;
		Peak	tPeak = new Peak(0f,0f);
			
		System.out.println("Normalizing spectrum...");
		//find max intensity in binnedSpectrum
		max = list.getMaxIntensity();
		//divide all bins by max and multiply by 100 (similar to XCorr)
		for (looper = 0; looper < binnedSpectrum.size(); looper++) {
                        normInt = (binnedSpectrum.get(looper).getIntensity()/max) * 100;
			tPeak = new Peak(binnedSpectrum.get(looper).getM2z(),normInt);
			binnedSpectrum.set(looper, tPeak);
                }
		return binnedSpectrum; 
	}

	 List<Peak> filterPeaks(ParamsFile params, PeakList list) {
                float   max = 0;
                int looper = 0;
                float   normInt;
                float   mz;
		List<Peak>      peakList;
		List<Peak>	tempList = new ArrayList();
                Peak    tPeak = new Peak(0f,0f);

                //sort peaks by intensity
                peakList = list.getSortedPeaks(true);
                //remove all peaks that are below intensity rank threshold
                if (peakList.size() - 1 >= params.PeakThreshold) { 
			for(looper = peakList.size() - 1 - params.PeakThreshold; looper < peakList.size() - 1; looper++) {
                        	tempList.add(peakList.get(looper));
                	}
		}
                return tempList;
        }
}
