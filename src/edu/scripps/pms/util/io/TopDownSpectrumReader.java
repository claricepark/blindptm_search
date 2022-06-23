/**
 * @file TopDownSpectrumReader.java
 * This is the source file for edu.scripps.pms.util.spectrum.TopDownSpectrumReader
 * @author Tao Xu
 * @date $Date
 */



package edu.scripps.pms.util.io;

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.mspid.MassSpecConstants;

public class TopDownSpectrumReader {

    private ArrayList<String> files; 
    
    private void getFiles(String dir, String extension) {
        
         files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(extension)) {
                 files.add(dir + "/" + s);
                 //System.out.println(s);
             }
         }
          
         System.out.println("num files added: " + files.size());
    }
    public TopDownSpectrumReader (String dir, String extension) throws IOException {
        getFiles(dir, extension);
    }
    public int getNumSpectra() throws IOException {
        return files.size();
    }
    public static void main(String args[]) throws Exception {
        String line = null;
        TopDownSpectrumReader sr = new TopDownSpectrumReader(args[0], args[1]);

        Iterator<PeakList> it = sr.getSpectra();
        int counter = 0;
        int numPeaks = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            PeakList list = it.next();
//System.out.println("RetentionTime: " + list.getRetentionTime());
//          System.out.println(list.getLoscan() + "\t" + list.getHiscan() + "\t");

//	    System.out.println(list.getSpectrumWithoutHlines());
	    Peak p;
	    StringBuffer sb = new StringBuffer();
	//System.out.println("===>>" + list.getZlines());

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

            //System.out.println(zline.getM2z() + " " + zline.getChargeState());

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

        System.out.println("Finished");
    }

    public ArrayList<PeakList> getSpectraList() throws IOException {
        Iterator<PeakList> it = getSpectra();
        ArrayList<PeakList> spectraList = new ArrayList(40000);
        while(it.hasNext()) {
            spectraList.add(it.next());
        }
        return spectraList;
    }
    public Iterator <PeakList> getSpectra() throws IOException {
        return new Iterator() {
            Iterator<String> filesIt = files.iterator();
            public boolean hasNext() {
                return filesIt.hasNext();
            }

            public Object next() {

                PeakList peaks = null;
                try {
String s = filesIt.next();
System.out.print("Processing: " + s);
                    peaks = getPeakList(s);
System.out.println("\tnumPeaks: " + peaks.numPeaks());
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



    private PeakList getPeakList(String file) throws IOException {
        PeakList peakList = new PeakList();
        BufferedReader br = new BufferedReader(new FileReader(file), 4096);
        br.readLine();
        int scanNo = Integer.parseInt(br.readLine().trim().split("\t")[1]);
        peakList.setHiscan(scanNo);
        peakList.setLoscan(scanNo);
         
        String line = null;
        while((line = br.readLine()) != null && !line.startsWith("Xtract Result")); 
        while((line = br.readLine()) != null) {
            String [] arr = line.trim().split("\t");
            if(arr.length == 8) {
                float mass = Float.parseFloat(arr[0]) + MassSpecConstants.MASSH;
                int charge = Integer.parseInt(arr[1]);
                peakList.addPeak(new Peak(mass, 1, charge));
            }
        }
                    
        br.close();
        return peakList;
    }

}
