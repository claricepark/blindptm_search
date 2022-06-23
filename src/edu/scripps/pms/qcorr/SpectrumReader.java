/**
 * @file SpectrumReader.java
 * This is the source file for edu.scripps.pms.util.spectrum.SpectrumReader
 * @author Tao Xu
 * @date $Date: 2007/02/27 01:31:02 $
 */



package edu.scripps.pms.qcorr;

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;

public class SpectrumReader {

    public static final char FIRSTCHAROFHLINE = 'H';

    // delimiter of m2z and intensity
    public static final String MZINTENSITYDELIMITER = " ";
    private BufferedReader br; // reads from the ms file
    // the format of ms data file, e.g., ms1, ms2, dta
    private String fileFormat;

    // if the spectraDelimiter is "S", then isNewFormat is true,
    // if the spectraDelimiter is ":", then isNewFormat is false
    private boolean isNewFormat = true;
    // delimiter for spectra, could be ':' or 'S', depends on the ms file type
    private char spectraDelimiter = 'S';
    private List<String> hlines = null;
    private String lastLine = ""; // remember the last line read

    public SpectrumReader (String msFileName, String fileFormat) throws IOException {
        br = new BufferedReader(new FileReader(msFileName), 4096);
        this.fileFormat = fileFormat;
        readHlines();
        if (lastLine != null) {
            spectraDelimiter = lastLine.charAt(0);
//System.out.println("spectraDelimiter: " + spectraDelimiter);
            if (spectraDelimiter == 'S') {
                isNewFormat = true;
            } else if (spectraDelimiter == ':') {
                isNewFormat = false;
            } else {
                throw new FileFormatUnknownException();
            }
        }
    }

    public static void main(String args[]) throws Exception {
        String line = null;
        SpectrumReader sr = new SpectrumReader(args[0], args[1]);
	Hline h = new Hline(sr.getHlines());

        Iterator<PeakList> it = sr.getSpectra();
        int counter = 0;
        int numPeaks = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            PeakList list = it.next();
//           System.out.println(list.getLoscan() + "\t" + list.getHiscan() + "\t");

//	    System.out.println(list.getSpectrumWithoutHlines());
	    Peak p;
	    StringBuffer sb = new StringBuffer();
	System.out.println("===>>" + list.getZlines());

	    for(Iterator<Peak> itr=list.getPeaks(); itr.hasNext(); )
	    {
		    p = itr.next();
		    sb.append(p.getM2z());
		    sb.append("\t");
		    sb.append(p.getIntensity());
		    sb.append("\n");
	    }



        for(Iterator itr = list.getZlines(); itr.hasNext(); )
        {
           Zline  zline = (Zline)itr.next();

		System.out.println(zline.getM2z() + " " + zline.getChargeState());

        }


 
            //list.sort(sortByIntensity);
//            System.out.println("MaxIntens: " + list.getMaxIntensity());
//            System.out.println("MinIntens: " + list.getMinIntensity());
            //System.out.print("Min m2z: " + list.getMinM2z());
            //System.out.print("\tMax m2z: " + list.getMaxM2z());
            //System.out.print("\tRange: " + (list.getMaxM2z()-list.getMinM2z()));

            //System.out.println("\tNum peaks: " + list.numPeaks());
            numPeaks += list.numPeaks();
            counter++;
            //if (counter > 20)
            //break;
        }
        System.out.println("Total number of spectra processed: " + counter);
        System.out.println("Average number of peaks per list: " + numPeaks/counter);
        sr.closeDataFile();

        System.out.println("Finished");
    }

    public Iterator <String> getHlines() {
        return hlines.iterator();
    }

    public Iterator <PeakList> getSpectra() throws IOException {
        return new Iterator() {
            	public boolean hasNext() {
                	return lastLine != null;
            	}

            	public Object next() {

                	PeakList peaks = null;
                	try {
                    		peaks = getPeakList();
                	} catch (IOException e) {
                    		e.printStackTrace();
                	}
                	return peaks;
            	}

            	public void remove() {
                	throw new UnsupportedOperationException("Not supported");
            	}
        };
    }

    public void closeDataFile() throws IOException {
        br.close();
    }

    private void readHlines() throws IOException {

        boolean isFirstRead = true; // if file is read for the first time
        if (hlines == null) {
            hlines = new ArrayList<String>();
        } else {
            isFirstRead = false;
        }
        while ((lastLine = br.readLine()) != null) {
            if (lastLine.charAt(0) == FIRSTCHAROFHLINE) {
                if (isFirstRead) {
                    hlines.add(lastLine);
                }
            } else if (lastLine.charAt(0) == 'S' || lastLine.charAt(0) == ':') {
                return;
            } else {
                continue;            
            }
        }
    }

    private PeakList getPeakList() throws IOException {
        PeakList peakList = new PeakList();
        peakList.setHlines(hlines);
        getListSpecificInfo(peakList);
        getPeaks(peakList);
        return peakList;
    }

    private void getPeaks(PeakList list) throws IOException {
        String [] content = null;
        while(lastLine != null && lastLine.charAt(0) != spectraDelimiter) {
	/*try {
            // process one line for m2z and intensity
            content = lastLine.split(MZINTENSITYDELIMITER);
            // add a new Peak to the peaklist
	    if (content.length < 2) {
            	list.addPeak(new Peak(Float.parseFloat(content[0]),
                                   Float.parseFloat(content[1])));
	    } else {
		list.addPeak(new Peak(Float.parseFloat(content[0]),
                                   Float.parseFloat(content[1]), Integer.parseInt(content[2])));
	    }
	} catch (Exception e) {System.out.println(e.toString() + lastLine);}
            // check next line*/
            try {
                // process one line for m2z and intensity
                content = lastLine.split(MZINTENSITYDELIMITER);
                // add a new Peak to the peaklist
                float m2z = Float.parseFloat(content[0]);
                float intens = Float.parseFloat(content[1]);
                if(m2z != 0 || intens != 0) { // avoid to add peaks that have m2z or intensity as 0
                    list.addPeak(new Peak(m2z, intens));
                }
            } catch (Exception e) {
               System.out.println(e.toString() + lastLine);
            }
                // check next line
            lastLine = br.readLine();
        }
    }
    private void getListSpecificInfo(PeakList list) throws IOException {

        String [] content = null;
        Zline z = null;
        if (isNewFormat){
            // parse S line
	
            content = lastLine.split("\t");
            list.setLoscan(Integer.parseInt(content[1]));
            list.setHiscan(Integer.parseInt(content[2]));
            if (fileFormat.equals("ms2")) {
                list.setPrecursorMass(Float.parseFloat(content[3]));
            }
            // parse Z lines
            while ((lastLine = br.readLine()) != null) {
                if (lastLine.charAt(0) == 'Z') {
                    content = lastLine.split("\t");
                    z = new Zline(Integer.parseInt(content[1]),
                                         Float.parseFloat(content[2]));
                    list.addZline(z);
                } else if (lastLine.charAt(0) == 'D') {
                    if(z != null) {
                        z.addDline(lastLine);
                    } else {
                        throw new RuntimeException("Wrong file format! Dline should always follow Zline");
                    }
               
                } else if (lastLine.charAt(0) == 'I') {
                    list.addIline(lastLine);
                } else if (lastLine.charAt(0) >= '0' && lastLine.charAt(0) <= '9') {
                    break; // digits, peaks starts
                } else {
                    continue;  // ignor lines other than I, D, S or peaks, e.g., empty line 
                }
            }
        } else {
            // process file in old format, i.e, spectra delimiter is ":"
           do {
                String line = lastLine.substring(1);
                content = line.split("\\.");

                list.setLoscan(Integer.parseInt(content[1]));
                list.setHiscan(Integer.parseInt(content[2]));
                if (fileFormat.equals("ms2")) {
  //                  list.setPrecursorMass(Float.parseFloat(content[3]));
                }
                 lastLine = br.readLine();
                 String [] oldZline = lastLine.split(" ");

                 z = new Zline(Integer.parseInt(oldZline[1]),
                                       Float.parseFloat(oldZline[0]));
                 list.addZline(z);
            } while ((lastLine = br.readLine()) != null && lastLine.charAt(0) == spectraDelimiter);
        }
    }
}
