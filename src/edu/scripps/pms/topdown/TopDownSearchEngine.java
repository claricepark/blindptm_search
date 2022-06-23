/**
 * @file TopDownSearchEngine.java
 * This is the source file for edu.scripps.pms.topdown.TopDownSearchEngine
 * @author Tao Xu
 * @date $Date: 2007/12/02 18:42:03 $
 */
package edu.scripps.pms.topdown;

import edu.scripps.pms.util.seq.Fasta;

import java.net.*;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.TopDownSpectrumReader;

class TopDownSearchEngine {

    //public static final int DEFAULTNUMPEPTIDEHIT = 1000;
    //public static final int DEFAULTPEPTIDELENGTH = 40;
    private SearchParams params;
    private ProteinDatabase pdb;
    private MassCalculator mc;
    //private ProcessedPeakList ppl = null;
    public static final int HYPERGEOMETRYSCORE = 1;
    private static final String USAGE = "java TopDownSearchEngine foldername";
    private String hostName;
    private long startTime;
    private long endTime;
    private double maxPrcMass;
    private double minPrcMass;

    private static StringBuffer HEADER = new StringBuffer(1200);
      /*  "H\tSQTGeneratorVersion\t 0.1\n"+
        "H\tComment\tThe Department of Cell Biology\n"+
        "H\tComment\tThe Scripps Research Institute\n" +
        "H\tComment\tProbability Score and Confidence developed\n" +
        "H\tComment\tby Rovshan Sadygov and JR Yates, III\n" +
        "H\tComment\tRef. Anal. Chem., 2003\n"+
        "H\tComment\tCan use Cross-Correlation of J.Eng and J.Yates\n" +
        "H\tComment Ref. J. Am. Soc. Mass Spectrom., 1994, 4, p. 976\n"+
        "H\tComment Paralellization program is prolucid_submit\n";
      */
    private static DistributionCalculator dc;

    public TopDownSearchEngine(SearchParams sp, ProteinDatabase db) { 
        this.params = sp;
        minPrcMass = params.getMinPrecursorMass(); 
        maxPrcMass = params.getMaxPrecursorMass();
        pdb = db;
        mc = new MassCalculator(sp);
        hostName = getHostName();
        getHeader();
    }

    private void getHeader() {
        HEADER.append("H\tProteinIDGenerator\tTopDownProLuCID\n");
        HEADER.append("H\tProteinIDGeneratorVersion\t0.1\n");
        HEADER.append("H\tDatabase\t" + params.getDbName()+ "\n");
        HEADER.append("H\tPrecursorMasses\t" + params.getPrecursorIsotope() + "\n");
        HEADER.append("H\tFragmentMasses\t" + params.getFragmentIsotope() + "\n");
    }
    public static String getHostName() {
        String hostName = null;
        try {
            hostName  = InetAddress.getLocalHost().getHostName();
            System.out.println("host name " + hostName);
            if(hostName.indexOf(".") != -1) {
                hostName = hostName.substring(0, hostName.indexOf("."));
            }
        } catch (UnknownHostException e) {
            e.printStackTrace();
        }
        return hostName;

    }
    public List<SearchResult> search(PeakList peaks) throws IOException {
        ArrayList<SearchResult> results = new ArrayList<SearchResult>();
        //PeptideHit.setFragMasses(mc.getFragMasses()); makes no big difference
            startTime = System.currentTimeMillis();
                double maxFragment = peaks.getMaxM2z() + 100; 
                ProcessedPeakList ppl = new ProcessedPeakList(peaks, new Zline(1, maxFragment), params, mc);
                SearchResult sr = pdb.topDownSearch(ppl);//search(peaks, params, z);
         
            //System.out.println("TimeUsed: " + timeUsedInSeconds);
                sr.calcScores();
                sr.setHostName(hostName);
            //System.out.println(sr.outputResults());
                //endTime = System.currentTimeMillis();
                //long timeUsedInSeconds = (endTime - startTime);
                sr.setSearchTime(System.currentTimeMillis() - startTime);
                results.add(sr);
            //System.out.println("ChargeState: " + chargeState + "\tPrecursorM2z: "
            //       + precMass + "\tNumber of peptideHits: " + sr.getNumPeptideHit());
        return results;
    }


    public static void main(String args[]) throws Exception {
        try {
            String dir = args[0];
            SearchParams sp = new SearchParams("topdown.xml");
            
            ProteinDatabase pd = new ProteinDatabase(sp.getDbName(), new MassCalculator(sp), sp);
            dc = new DistributionCalculator();
            //ProteinDatabase pd = new IndexedProteinDatabase(args[0]);
            TopDownSearchEngine se = new TopDownSearchEngine(sp, pd);
            TopDownSpectrumReader sr = new TopDownSpectrumReader(dir, "xls");
            String outFileName =  "proteinid.txt";
            PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName)));
 
            outFile.print(HEADER);
            //int numSearch = Integer.parseInt(args[2]);
            //System.out.println();
            int numPeaksMin = sp.getMinNumSpectra(); 
            int numPeaksMax = sp.getMaxNumSpectra(); 
            long beginTime = System.currentTimeMillis();
            int i = 0;
            TimeUtils timer = new TimeUtils();
            timer.startTiming();
            for (Iterator<PeakList> it = sr.getSpectraList().iterator(); it.hasNext();) {
                i++;
                PeakList pl = it.next();        
                int numPeaks = pl.numPeaks();
                if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax) {
                //timer.stopTiming(); 
                //System.out.println("Time used to read spectrum number " + i + ": " + timer.getTimeUsedMillis());
                //double prcMass = pl.getPrecursorMass();
                    List<SearchResult> results = se.search(pl);
                    for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {
                        
                        outFile.print(resultIt.next().outputResults());
                    //outFile.flush();
                    }
               
                    System.gc();
                } 


                //timer.startTiming();
                //if (i == numSearch) {
                //    break;
               // }
            }
            long endTime = System.currentTimeMillis();
            long timeUsed = endTime - beginTime;
            outFile.close();
            System.out.println("Time used for search: " + timeUsed);
            System.out.println();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
}


