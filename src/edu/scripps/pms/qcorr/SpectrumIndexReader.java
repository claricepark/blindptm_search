/**
 * @file SpectrumIndexReader.java
 * @author Robin Park
 * @date $Id: SpectrumIndexReader.java,v 1.3 2008/05/07 22:36:27 rpark Exp $
 */



package edu.scripps.pms.qcorr;

import java.io.*;
import java.util.*;

public class SpectrumIndexReader {

//    public static final char FIRSTCHAROFHLINE = 'H';
    private RandomAccessFile file=null;

    // delimiter of m2z and intensity
//    public static final String MZINTENSITYDELIMITER = " ";
    private BufferedReader br; // reads from the ms file
    // the format of ms data file, e.g., ms1, ms2, dta
    private String spectrumNum;
    // if the spectraDelimiter is "S", then isNewFormat is true,
    // if the spectraDelimiter is ":", then isNewFormat is false
    //private boolean isNewFormat = true;
    // delimiter for spectra, could be ':' or 'S', depends on the ms file type
    private final String SLINE_HEADER = "S";
    private final String ZLINE_HEADER = "Z";
    //private List<String> hlines = null;
    //private String lastLine = ""; // remember the last line read

    public SpectrumIndexReader (String fileName, String spectrumNum) throws IOException {

	this.spectrumNum = spectrumNum;
	try
	{
		br = new BufferedReader(new FileReader(fileName + ".index"), 4096);
	}
	catch(IOException e)
	{
		System.out.println("index file does not exist " + e);
		throw new IOException("index file does not exist " + e);
	}

       	file = new RandomAccessFile(fileName, "r");
    }

    public PeakList getPeakList() throws IOException
    {
        String[] strArr = getPeaks().split("\n");
        
        Peak peak;
        PeakList peakList = new PeakList();
        for(int i=1;i<strArr.length;i++)
        {
            String[] peakArr = strArr[i].split(" ");

            peakList.addPeak( new Peak(Float.parseFloat(peakArr[0]), Float.parseFloat(peakArr[1])) );

        }

        return peakList;
    }

    public PointList getPointList() throws IOException
    {
        String[] strArr = getPeaks().split("\n");

        Point point;
        PointList pointList = new PointList();
        for(int i=1;i<strArr.length;i++)
        {
            String[] pointArr = strArr[i].split(" ");

            pointList.addPoint( new Point(Float.parseFloat(pointArr[0]), Float.parseFloat(pointArr[1])) );
        }
        return pointList;
    }

    public String getClosestScan(long scan) throws IOException
    {
	String lastLine;
        String st=null;

	while ((lastLine = br.readLine()) != null) {
		lastLine = lastLine.substring(0, lastLine.indexOf("\t"));
		System.out.println("lastLine == " + lastLine);
		if (scan <= Long.parseLong(lastLine)) {
			st = lastLine;
		}
	}

        if(null == lastLine)
            throw new IOException("Couldn't locate closest spectrum : ");

	return st;
    }

    public String getPeaks() throws IOException
    {
        String lastLine;
        StringBuffer sb = new StringBuffer();
	//System.out.println("Looking for spectrum in index file...");
        while ((lastLine = br.readLine()) != null && !lastLine.startsWith(spectrumNum));

        if(null == lastLine)
            throw new IOException("Spectrum number does not exist : ");

        lastLine = lastLine.substring(lastLine.indexOf("\t") +1);

        file.seek(Long.parseLong(lastLine));

        //Skip S line
        file.readLine();

	//Skip Z lines
	lastLine = file.readLine();
	while(lastLine.startsWith(ZLINE_HEADER))
	{
		lastLine = file.readLine();
	}

	//Read Spectra file
	while(lastLine != null && !lastLine.startsWith(SLINE_HEADER))
	{
		sb.append(lastLine);
		sb.append("\n");

		lastLine = file.readLine();
	}

	return sb.toString();
    }

    public void close() throws IOException {
        file.close();
        br.close();
    }


    public static void main(String args[]) throws Exception
    {

        SpectrumIndexReader reader = new SpectrumIndexReader(args[0], args[1]);

        for(ListIterator<Peak> itr=reader.getPeakList().getPeaks(); itr.hasNext(); )
        {
            Peak p = itr.next();

            System.out.println(p + "\t" + p.getM2z() + "\t" + p.getIntensity());
        }
    }
       
}
