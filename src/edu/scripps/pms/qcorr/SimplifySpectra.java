package edu.scripps.pms.qcorr;

import java.util.*;
import java.io.*;
import java.text.*;

public class SimplifySpectra
{
	static float           carbon13Mass = 1.003f;
	static float           massTolerance = 100f;//units of ppm
        static int             maxChg = 8;

	public PeakList simplifySpectra(PeakList list, float massTolerance, boolean profileMode) {
		PeakList        pickedList = list;
		int numpeaks = 0;

		SimplifySpectra simplify = new SimplifySpectra();
		System.out.println("Picking Peaks....");
		if (profileMode) {
			pickedList = simplify.peakPicker(list);
		}
		numpeaks = pickedList.numPeaks();
		//simplify.printLists(list, pickedList);
		if (numpeaks > 3) {
			System.out.println("Assigning Charge State....");
			simplify.assignChargeState(pickedList);
		//	simplify.printLists(list, pickedList);
			System.out.println("Reassigning MZ....");
			pickedList = simplify.reassignMZ(pickedList);
		}
		simplify.printLists(list, pickedList);
		return pickedList;
	}

	public static void main(String[] args) throws Exception {
	int             dirIterator;
	String		currentFile;
	PeakList 	list;
	Peak 		p;
	int		scan;
	float		precursor;
	boolean		makeMS2 = true;
	String		temp;
	boolean		profileMode = false;

	//for(dirIterator = 0; dirIterator < args.length; dirIterator++) {
        	//currentFile = args[dirIterator];
		currentFile = args[0];
		System.out.println(args[0] + "\t" + args[1] + "\t" + args[2]);
		massTolerance = new Float(args[1]).floatValue();
		maxChg = new Integer(args[2]).intValue();
                System.out.println("Reading " + currentFile + "...");
                SpectrumReader ms2 = new SpectrumReader(currentFile, "ms2");
		Iterator<PeakList> spectrumList = ms2.getSpectra();
		Iterator<String> hLines = ms2.getHlines();
		int counter = 0;
		int numPeaks = 0;
		int numpeaks;
		SimplifySpectra      simplify = new SimplifySpectra();
		if (makeMS2) {
			simplify.printHLines(currentFile, hLines);
		}
		while (spectrumList.hasNext()) {
			list = spectrumList.next();
			PeakList	pickedList = list;
			System.out.println("Processing Spectrum" + "\t" + list.getLoscan());
			System.out.println("Picking Peaks....");
			if (profileMode) {
                        	pickedList = simplify.peakPicker(list);
                	}
			//pickedList = simplify.peakPicker(list);
			numpeaks = pickedList.numPeaks();
	                if (numpeaks > 3) {
				if (makeMS2) {
                                	simplify.printSandZLines(currentFile, list);
                        	}
				System.out.println("Assigning Charge State....");
				simplify.assignChargeState(pickedList);
				System.out.println("Reassigning MZ....");
				pickedList = simplify.reassignMZ(pickedList);
				if (makeMS2) {
					Iterator mzlist = pickedList.getPeaks();
					while (mzlist.hasNext()) {
						p = (Peak)mzlist.next();
						simplify.writeOutputLine("Simplified_" + currentFile, p.getM2z() + " " + p.getIntensity() + "\n", true);
					}
				}
				simplify.printSpectrum(pickedList, true, list.getLoscan());
		    		numPeaks += pickedList.numPeaks();
		    		counter++;
			}
		}
		System.out.println("Total number of spectra processed: " + counter);
		System.out.println("Average number of peaks per list: " + numPeaks/counter);
		ms2.closeDataFile();
		System.out.println("Finished");
	//}
	}

	public void printHLines (String currentFile, Iterator<String> hLines) {
		String 	temp;
		if (hLines.hasNext()) {
			temp = hLines.next() + "\n";
			writeOutputLine("Simplified_" + currentFile, temp, false);
		}
		while(hLines.hasNext()) {
                        temp = hLines.next() + "\n";
                	writeOutputLine("Simplified_" + currentFile, temp, true);
                }
	}

	public void printSandZLines (String currentFile, PeakList list) {
                String  temp;
        	String  templowscan = "000000";
        	String  temphighscan = "000000";
		String	sLine;
		String	zLine;
		int	chgState;
		float	mzMass;

		if (list.getLoscan() < 10) {
			templowscan = "00000" + list.getLoscan();
			temphighscan = "00000" + list.getHiscan();
		} else if (list.getLoscan() >= 10 && list.getLoscan() < 100) {
			templowscan = "0000" + list.getLoscan();
			temphighscan = "0000" + list.getHiscan();
		} else if (list.getLoscan() >= 100 && list.getLoscan() < 1000) {
			templowscan = "000" + list.getLoscan();
			temphighscan = "000" + list.getHiscan();
		} else if (list.getLoscan() >= 1000 && list.getLoscan() < 10000) {
			templowscan = "00" + list.getLoscan();
			temphighscan = "00" + list.getHiscan();
		} else if (list.getLoscan() >= 1000 && list.getLoscan() < 10000) {
			templowscan = "0" + list.getLoscan();
			temphighscan = "0" + list.getHiscan();
		} else {
			templowscan = Integer.toString(list.getLoscan());
			temphighscan = Integer.toString(list.getHiscan());
		}

		sLine = "S" + "\t" + templowscan + "\t" + temphighscan + "\t" + list.getPrecursorMass() + "\n";
		writeOutputLine("Simplified_" + currentFile, sLine, true);
		Iterator<Zline> zLines = list.getZlines();
		while (zLines.hasNext()) {
			Zline temp2 = zLines.next();
			chgState = temp2.getChargeState();
			mzMass = temp2.getM2z();
			zLine = "Z" + "\t" + chgState + "\t" + mzMass + "\n";
			writeOutputLine("Simplified_" + currentFile, zLine, true);
		}
        }
	
	public PeakList peakPicker(PeakList list) {
		Peak            p;
		Peak		pAfter;
		Peak		pBefore;
		Iterator 	peakList;
		PeakList	pickedList = new PeakList();
		int		numpeaks;
        
		peakList = list.getPeaks();
		numpeaks = list.numPeaks();
		if (numpeaks > 3) {
			pBefore = (Peak)peakList.next();
			p = (Peak)peakList.next();
			for(list.getPeaks(); peakList.hasNext(); ) {
				pAfter = (Peak)peakList.next();
				if (p.getIntensity() > pAfter.getIntensity() && p.getIntensity() > pBefore.getIntensity()) {
                                	pickedList.addPeak(p);
				}
				pBefore = p;
				p = pAfter;
			}
			copyVariables(pickedList, list);
			return pickedList;
		}
		else {
			return list;
		}
	}

	public void copyVariables(PeakList pickedList, PeakList list) {
		pickedList.setLoscan(list.getLoscan());
		pickedList.setHiscan(list.getHiscan());
		Iterator<Zline> temp = list.getZlines();
		ArrayList<Zline> temp2 = new ArrayList();
		while(temp.hasNext()) {
			temp2.add(temp.next());
		}
		pickedList.setZlines(temp2);
		pickedList.setPrecursorMass(list.getPrecursorMass());
		pickedList.setListType(list.getListType());
	}

	public float getChargeDifference(float mz1, float mz2) {
	    int charge = maxChg;
	    int testValue; 
	    float floatCharge;

	    while (charge > 0) {
		floatCharge = new Float(charge).floatValue();
		testValue = this.getMassDifference(mz1, mz2, charge);
		if (testValue == -1)
		    return (floatCharge + 0.5f);
		else if (testValue == 0)
		    return floatCharge;
		else
		    charge--;
	    }
	    return 0.5f;
	}

        public int getMassDifference(float mz1, float mz2, int charge) {
            float floatCharge = new Float(charge).floatValue();
	    float deltaMass = carbon13Mass/floatCharge;
	    float massDifference;
	    float upperMass;
	    float lowerMass;

	    if (mz1 < mz2)
		massDifference = mz2 - mz1;
	    else
		massDifference = mz1 - mz2;

	    upperMass = massDifference + (mz1+mz2)*massTolerance/2000000.0f;
	    lowerMass = massDifference - (mz1+mz2)*massTolerance/2000000.0f;
	    if (deltaMass > upperMass)
		return -1;
	    else if (deltaMass < lowerMass)
		return 1;
	    else
		return 0;
        }

	public boolean isIntegerCharge(float testCharge) {
	    int charge = maxChg;
	    float floatCharge;

	    while (charge > 0) {
                floatCharge = new Float(charge).floatValue();
                if (testCharge > floatCharge)
                    return false;
                else if (testCharge == floatCharge)
                    return true;
                else
                    charge--;
            }
            return false;

	}

	public void assignChargeState(PeakList list) {
		Peak            peakRunner;
		Peak		peak2;
		ListIterator 	peakList;
		float		mz;
		float		floatCharge;
		float		currentCharge;
		int		charge;
		boolean		stopCondition;
		int		peakCounter;
		int		tempCounter;
		int 		numPeaks;

		peakList=list.getPeaks();
		numPeaks  = list.numPeaks();
                if (numPeaks > 3) {
			peakRunner = (Peak)peakList.next();
			while (peakList.hasNext()) {
			    if (peakRunner.getChargeState() == 0) {
				mz = peakRunner.getM2z();
				stopCondition = false;
				peak2 = peakRunner;
				peakCounter = 0;
				while (peakList.hasNext() && stopCondition == false) {
				    peak2 = (Peak) peakList.next();
				    peakCounter++;
				    if (peak2.getChargeState() == 0) {
					floatCharge = this.getChargeDifference(mz, peak2.getM2z());
					if (peakRunner.getChargeState() == 0) {
					    if (this.isIntegerCharge(floatCharge)) {
						charge = new Float(floatCharge).intValue();
						peakRunner.setChargeState(charge);
						peak2.setChargeState(charge);
						mz = peak2.getM2z();
					    }
					    else if (floatCharge == 0.5f) {
						stopCondition = true;	
					    }
					}
					else {
					    currentCharge = new Float(peakRunner.getChargeState()).floatValue();
					    if (floatCharge == currentCharge) {
						charge = new Float(floatCharge).intValue();
						peak2.setChargeState(charge);
						mz = peak2.getM2z();
					    }
					    else {
						if (floatCharge < currentCharge) {
						    stopCondition = true;
						}
						if (this.isIntegerCharge(floatCharge)) {
						}
					    }
					}
				    }
				}
				for (tempCounter = 0; tempCounter < peakCounter; tempCounter++) {
				    peak2 = (Peak) peakList.previous();
				}
			    }
			    peakRunner = (Peak) peakList.next();
			}
		}
	}


	public PeakList reassignMZ(PeakList pickedList) {
		Peak            p;
		int		chg;
		Iterator        peakList;
                PeakList        reAssignedList = new PeakList();

		peakList = pickedList.getPeaks();
		p = (Peak)peakList.next();
		for(pickedList.getPeaks(); peakList.hasNext(); ) {
                	p = (Peak)peakList.next();
                        if (p.getChargeState() > 1) {
				chg = p.getChargeState();
				p.setM2z((p.getM2z() * chg) - ((chg - 1) * 1.00728f));
				reAssignedList.addPeak(p);
                        }
			else {
				reAssignedList.addPeak(p);
			}
                }
		//sort by mz
		reAssignedList.sortPeaks(false);
		copyVariables(reAssignedList, pickedList);
		return reAssignedList;
	}

	public void printLists(PeakList list, PeakList pickedList) {
        	Peak            p;
		Iterator peakList;

		//print original data
		System.out.println("Original Data");
		for(peakList=list.getPeaks(); peakList.hasNext(); ) {
			//Print specs
			StringBuffer sb = new StringBuffer();
                        p = (Peak)peakList.next();
                        sb.append(p.getM2z());
                        sb.append("\t");
                        sb.append(p.getIntensity());
			sb.append("\t");
			sb.append(p.getChargeState());
                        System.out.println(sb);
                }
		//print peak picked and charge assigned data
		System.out.println("Simplified Data");
		for(peakList=pickedList.getPeaks(); peakList.hasNext(); ) {
                        //Print specs
                        StringBuffer sb = new StringBuffer();
                        p = (Peak)peakList.next();
                        sb.append(p.getM2z());
                        sb.append("\t");
                        sb.append(p.getIntensity());
			sb.append("\t");
                        sb.append(p.getChargeState());
                        System.out.println(sb);
                }
        }
	
	
	public void printSpectrum(PeakList spectrum, boolean OverWrite, int specnum) {
		String          Filename = "Simplify.txt";
		try {
                        File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);
                	Iterator        peakList;
			//peakList = spectrum.getPeaks();
			Peak		p;
			Outgoing.write("Spectrum\t" + specnum + "\n");
                	for(peakList = spectrum.getPeaks(); peakList.hasNext(); ) {
				String	Data = new String();
                        	p = (Peak)peakList.next();
				Data = p.getM2z() + " " + p.getIntensity() + "\n";// + "\t" + p.getChargeState() + "\n";
				Outgoing.write(Data);
                        }
			Outgoing.close();
                } catch (IOException failure) {
                        System.out.println("IO Error while writing " + Filename);
                        System.exit(0);
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
}
