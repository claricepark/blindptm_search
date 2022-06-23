/**
 * @file ProlucidSearchEngine.java
 * This is the source file for edu.scripps.pms.mspid.ProlucidSearchEngine
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.mspid;

import java.net.*;
import java.util.*;
import java.io.*;
import java.lang.management.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.mspid.db.MemorizedIndexedProteinDatabase;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.MzxmlSpectrumReader;
import org.jdom.JDOMException;


// mutli threading prolucid
public class ProlucidSearchEngine {

    //public static final int DEFAULTNUMPEPTIDEHIT = 1000;
    //public static final int DEFAULTPEPTIDELENGTH = 40;
    private SearchParams params;
    private ProteinDatabase pdb;
    private String spectrumFile;
    private MassCalculator mc;
    //private ProcessedPeakList ppl = null;
    public static final int HYPERGEOMETRYSCORE = 1;
    private static final String USAGE = "Usage: java -Xmx500M ProlucidSearchEngine spectrumfilenameOrfolder searchparameterfile max_num_thread";
    private String hostName;
    //private long startTime;
    private long endTime;
    private double maxPrcMass;
    private double minPrcMass;

    //private static StringBuffer HEADER = new StringBuffer(1200);
    //private static StringBuffer HEADER = null;
    //private StringBuffer HEADER = null;
    private static DistributionCalculator dc;
    private Iterator<PeakList> itr = null;
    private int numSpectraSearched = 0;
    private int numSpectraIgnored = 0;
    private PrintWriter outFile = null;
    private int numPeaksMin = 0; 
    private int numPeaksMax = 0; 
    private int numSpectraInFile = 0;
    protected static StringBuffer logbuffer = new StringBuffer(2000);
    private int ignoredResults =0;
    private int nResults =0;
    private int nMinPrcMass =0;
    private int nPtrue =0;
    private int nMaxPrcMass =0;


    public ProlucidSearchEngine(SearchParams sp, ProteinDatabase db, String msFileName) 
                                                   throws IOException, JDOMException, Exception{ 
        this.params = sp;
        spectrumFile = msFileName;
        minPrcMass = params.getMinPrecursorMass(); 
        maxPrcMass = params.getMaxPrecursorMass();
        pdb = db;
        mc = new MassCalculator(sp);
        hostName = getHostName();
        String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".")) + ".sqt";
        // default buffer size is 8192, 40960 is 5X of the default size
        outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName), 40960));
        //if(HEADER == null) {
        //    getHeader();
       // }
        //outFile.print(HEADER);
        outFile.print(getHeader());
        ArrayList<PeakList> spectra = getSpectra(msFileName);
        numSpectraInFile =  spectra.size();
        itr = spectra.iterator();
        numPeaksMin = sp.getMinNumSpectra(); 
        numPeaksMax = sp.getMaxNumSpectra(); 
        if(dc == null) {
            dc = new DistributionCalculator();
        }

    }

    public synchronized PeakList getSpectrum() {
        if(itr != null && itr.hasNext()) {
            return itr.next(); 
        }
        return null; 
    }

    private String getHeader() {
        StringBuffer HEADER = new StringBuffer(1200);
        HEADER.append("H\tSQTGenerator\tProLuCID\n");
        HEADER.append("H\tSQTGeneratorVersion\t1.4\n");
        HEADER.append("H\tComment ProlucidVersion\t1.8\n");
        HEADER.append("H\tComment ProLuCID is developed in the Yates laboratory at The Scripps Research Institute, La Jolla, CA\n");
        HEADER.append("H\tComment ProLuCID ref. Xu T, Venable JD, Park SK, Cociorva D, Lu B, Liao L, Wohlschlegel J, Hewel J, Yates JR 3rd\n");
        HEADER.append("H\tComment ProLuCID ref. ProLuCID, a fast and sensitive tandem mass spectra-based protein identification program.\n");
        HEADER.append("H\tComment ProLuCID ref. MOL CELL PROTEOMICS vol. 5(10): S174-S174 671 Suppl. S OCT 2006\n");
        HEADER.append("H\tComment Paralellization Program using PBS is submit_prolucid\n");
        HEADER.append("H\tComment Please send bug report or comments to Tao Xu by email taoxu@scripps.edu\n");
        HEADER.append("H\tHostName\t" + getHostName() + "\n");
        HEADER.append("H\tSpectrumFile\t" + spectrumFile + "\n");
        HEADER.append("H\tDatabase\t" + params.getDbName()+ "\n");
        HEADER.append("H\tDBLocusCount\t" + pdb.getNumSequences() + "\n");
        HEADER.append("H\tNumOutput\t" + params.getNumOutput() + "\n");
        HEADER.append("H\tPrecursorMasses\t" + params.getPrecursorIsotope() + "\n");
        HEADER.append("H\tFragmentMasses\t" + params.getFragmentIsotope() + "\n");
        int numIsotopicPeaks = params.getNumIsotopicPeaks();
        if(numIsotopicPeaks == 0) {
            HEADER.append("H\tHighPrecursorMassTolerance\t" + params.getHighPrecursorTolerance() + "\n");
            HEADER.append("H\tLowPrecursorMassTolerance\t" + params.getLowPrecursorTolerance() + "\n");
        } else {
            HEADER.append("H\tNumPrecursorIsotopicPeaks\t" + numIsotopicPeaks + "\n");
            HEADER.append("H\tPrecursorMassTolerance\t" + params.getPrecursorTolerance() + "\n");
        }
        HEADER.append("H\tFragmentMassTolerance\t" + params.getFragmentTolerance() + "\n");
        if(params.getStaticNTermMod() != 0) {
            HEADER.append("H\tNTermStaticMod\t" + params.getStaticNTermMod() + "\n");
        }
        if(params.getStaticCTermMod() != 0) {
            HEADER.append("H\tCTermStaticMod\t" + params.getStaticCTermMod() + "\n");
        }
        for(Iterator<Modification> it = params.getStaticMods(); it.hasNext();) {
            Modification m = it.next();
            HEADER.append("H\tStaticMod\t" + m.getResidue()+ "=" + mc.getPrecursorMass(m.getResidue()) + "\n");
        }
        
        if(params.getNumNTermDiffMods() > 0) {
            HEADER.append("H\tNTermDiffMod");
            for(Iterator<TerminalModification> it =  params.getNTermDiffMods(); it.hasNext();) {
                TerminalModification m = it.next();
                HEADER.append("\t" + m.getMassShift());
            }
            HEADER.append("\n");
        }
 
        if(params.getNumCTermDiffMods() > 0) {
            HEADER.append("H\tCTermDiffMod");
            for(Iterator<TerminalModification> it =  params.getCTermDiffMods(); it.hasNext();) {
                TerminalModification m = it.next();
                HEADER.append("\t" + m.getMassShift());
            }
            HEADER.append("\n");
        }

        HEADER.append("H\tMaxNumInternalDiffModsPerPeptide\t" + params.getMaxAlter() + "\n");
        if(params.getMaxAlter() > 0) {
            for(Iterator<DiffMod> it = params.getDiffMods(); it.hasNext();) {
                DiffMod m = it.next();
                //HEADER.append("H\tInternalDiffMod\t" + m.toString() + "\n");
                HEADER.append("H\tDiffMod\t" + m.toString() + "\n");
            }
        }
        int enzymespec = params.getEnzymeSpecificity();
        if(enzymespec > 0) {
            HEADER.append("H\tEnzymeSpecificity\t" + params.getEnzymeSpecificity() + "\n");
            Protease p = params.getEnzyme();
            boolean isC = p.getType();

            HEADER.append("H\tEnzymeName\t" + p.getName() + "\n");
            if(isC) {
                HEADER.append("H\tEnzymeEnd\tCTerm\n");
            } else {
                HEADER.append("H\tEnzymeEnd\tNTerm\n");
            }
            HEADER.append("H\tEnzymeResidues\t");
            for(int i = 0; i < MassSpecConstants.NUMCHARS; i++) {
                if(p.isDigestable((char)i)) {
                    HEADER.append((char)i); 
                }
            }
            HEADER.append("\n");

        } else {
            HEADER.append("H\tEnzymeSpecificity\tNo_Enzyme\n");
        }
        return HEADER.toString();
    }
    public static String getHostName() {
        String hostName = null;
        try {
            hostName  = InetAddress.getLocalHost().getHostName();
            //logbuffer.append("host name " + hostName + "\n");
            if(hostName.indexOf(".") != -1) {
                hostName = hostName.substring(0, hostName.indexOf("."));
            }
        } catch (UnknownHostException e) {
            //e.printStackTrace();
        }
        if(hostName == null) {
            hostName =  "UnknowHost";
        }
        return hostName;
    }


    public  ArrayList<SearchResult>  search(PeakList peaks) throws IOException {
        ArrayList<SearchResult> results = new ArrayList<SearchResult>();
        int numPeaks = peaks.numPeaks();
//System.out.println("Now processing scan " + peaks.getLoscan());
        StringBuffer mylogline = new StringBuffer();

        if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax) {
            if("HCD".equals(peaks.getActivationType())) {
                //peaks.processHcdSpectrum();
                //System.out.println("With HCD processing, number of peaks before: " + numPeaks + " and after: " +  peaks.numPeaks());
            }
            //System.out.print("NumSearched\t" + (numSpectraSearched+numSpectraIgnored) + "\t" + numSpectraInFile + "\t" );
            //PeptideHit.setFragMasses(mc.getFragMasses()); makes no big difference

            //String logline = "NumSearched\t" + totalNumProcessedSpectra + "\t" + numSpectraInFile + "\t" ;
            //String logline = "NumSearched\t" + totalNumProcessedSpectra + "\t" + numSpectraInFile + "\tScan#: " + peaks.getLoscan() + "\t" ;
            //logbuffer.append(logline);
            for(Iterator<Zline> zlines = peaks.getZlines(); zlines.hasNext();) {
                Zline z = zlines.next();
                double precMass = z.getM2z();
                int chargeState = z.getChargeState(); 
                //int numPeaks = peaks.numPeaks();
                if(chargeState < params.getMinPrecursorCharge() || chargeState > params.getMaxPrecursorCharge()) {
                    SearchResult sr = new SearchResult(z,peaks);
                    long  startTime = System.nanoTime();
                    sr.setHostName(hostName);
                    sr.setPrintOnlySLine(true);
                    //System.out.println(sr.outputResults());
                    //endTime = System.nanoTime();
                    //long timeUsedInSeconds = (endTime - startTime);
                    sr.setSearchTime((System.nanoTime() - startTime)/1000000);
                    results.add(sr);
                    nResults++;
                    continue;
                }
        
                ProcessedPeakList ppl = new ProcessedPeakList(peaks, z, params, mc); 
                //logbuffer.append("FNumPeaks: " + ppl.getFinalNumPeaks() + "\tScan#: " + peaks.getLoscan() + "\t+" + chargeState + "\tNumPeaks: " + peaks.numPeaks() + "\tM+H+: " + precMass + "\tptrue: " + ppl.getPTrue() + "\n");
                mylogline.append("\tFNumPeaks: " + ppl.getFinalNumPeaks() + "\t+" + chargeState + "\tNumPeaks: " + peaks.numPeaks() + "\tM+H+: " + precMass + "\tptrue: " + ppl.getPTrue());

                if(precMass >= minPrcMass && precMass <= maxPrcMass && ppl.getPTrue() < 0.5) {
                    long  startTime = System.nanoTime();
                    //SearchResult sr = ppl.search(pdb);//search(peaks, params, z);
                    SearchResult sr = pdb.search(ppl);//search(peaks, params, z);
            
                //System.out.println("TimeUsed: " + timeUsedInSeconds);
                    sr.calcScores();
                    sr.setHostName(hostName);
                //System.out.println(sr.outputResults());
                    //endTime = System.nanoTime();
                    //long timeUsedInSeconds = (endTime - startTime);
                    sr.setSearchTime((System.nanoTime() - startTime)/1000000);
                    results.add(sr);
                    nResults++;
                }
                else
                {
                    //if(precMass < minPrcMass) nMinPrcMass++;
                    //if(precMass > maxPrcMass) nMaxPrcMass++;
                    //if(ppl.getPTrue()>=0.5) nPtrue++;
                    //ignoredResults++;
                    SearchResult sr = new SearchResult(ppl);
                    long  startTime = System.nanoTime();
                    sr.setHostName(hostName);
                    sr.setPrintOnlySLine(true);
                    //System.out.println(sr.outputResults());
                    //endTime = System.nanoTime();
                    //long timeUsedInSeconds = (endTime - startTime);
                    sr.setSearchTime((System.nanoTime() - startTime)/1000000);
                    results.add(sr);
                    nResults++;
                }
                //System.out.println("ChargeState: " + chargeState + "\tPrecursorM2z: "
                //       + precMass + "\tNumber of peptideHits: " + sr.getNumPeptideHit());
            }
            if(params.getChargeDisambiguation()) {
                results = chargeDisambiguation(results);
            }
            
            //outputSearchResults(results);
            numSpectraSearched++;
   
           // System.gc();
            //Thread.currentThread ().yield ();        

        } else {
            numSpectraIgnored++;
            for(Iterator<Zline> zlines = peaks.getZlines(); zlines.hasNext();) {
                Zline z = zlines.next();
                SearchResult sr = new SearchResult(z,peaks);
                long  startTime = System.nanoTime();
                sr.setHostName(hostName);
                sr.setPrintOnlySLine(true);
                //System.out.println(sr.outputResults());
                //endTime = System.nanoTime();
                //long timeUsedInSeconds = (endTime - startTime);
                sr.setSearchTime((System.nanoTime() - startTime)/1000000);
                results.add(sr);
                nResults++;
            }
        } 

        //mylogline.append("\n");
        mylogline.append("\t" + (Runtime.getRuntime().totalMemory()/1000000 + "\tFree memory\t" + Runtime.getRuntime().freeMemory()/1000000) + "\n");
        String logline = "NumSearched\t" + (numSpectraSearched + numSpectraIgnored) + "\t" + numSpectraInFile + "\tScan#: " + peaks.getLoscan() ;
        outputLog(logline + mylogline.toString());


        //System.out.println();
        return results;
    }

    
    public ArrayList<SearchResult> chargeDisambiguation(ArrayList<SearchResult> results) {

        ArrayList<SearchResult> newResults = new ArrayList<SearchResult>();
        SearchResult topPrimary = null;
        SearchResult topSp = null;

        for(Iterator<SearchResult> it = results.iterator(); it.hasNext();) {
            SearchResult sr = it.next();
            if(topSp == null) {
                topSp = sr;
                topPrimary = sr;
            } else {
                if(sr != null) {
                    if(topSp.getTopHit().getSecondaryScore() < sr.getTopHit().getSecondaryScore()) {
                        topSp = sr;
                    }
                    if(topPrimary.getTopHit().getPrimaryScore() < sr.getTopHit().getPrimaryScore()) {
                        topPrimary = sr;
                    }
                }
            }
        }
        if(topSp != null) {
            newResults.add(topPrimary);
            if(topSp != topPrimary) {
                newResults.add(topSp);
            } 
        }

        return newResults;

    }
    private static ArrayList<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception {
        ArrayList<PeakList> peaklists = null;
        if(file.endsWith(".mzXML")) {
            MzxmlSpectrumReader sr = new MzxmlSpectrumReader(file);
            ArrayList<PeakList> spectra = new ArrayList<PeakList>(20000);
            for(Iterator<MzxmlPeakList> it = sr.getSpectra(2); it.hasNext();) {
                spectra.add(it.next()); 
            }
            peaklists = spectra;
            sr.closeDataFile();
            sr = null;
        } else {
            SpectrumReader sr = new SpectrumReader(file, "ms2");
            peaklists = sr.getSpectraList();
            sr.closeDataFile();
            sr = null;
        }
    
        return peaklists;
    }
    public void closeOutputFile() throws IOException {

        outFile.close();
    }
    private synchronized void outputLog(String logline) {
        logbuffer.append(logline);
        int totalNumProcessedSpectra = numSpectraSearched+numSpectraIgnored; 
        if(totalNumProcessedSpectra%20 == 0) {
            System.out.print(logbuffer.toString());
            logbuffer = new StringBuffer(2000);
        } 
    }
    private synchronized void outputSearchResults(List<SearchResult> results) {
        
        for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {
            
            SearchResult rest = resultIt.next();
            if(rest != null) { 
                outFile.print(rest.outputResults());
            }
            //outFile.flush();
        }

    }
    public synchronized void outputSearchResults(StringBuffer results) {
        
            
        outFile.print(results.toString());
            //outFile.flush();

    }
    public static ArrayList<String> getMs2Files(String f) {
        ArrayList<String> ms2files = new ArrayList<String>();
        File myfile = new File(f);
        if(myfile.exists()) {
            if(myfile.isDirectory()) {
                String [] allfiles = myfile.list();
                for(int i = 0; i < allfiles.length; i++) {
                    String file = allfiles[i];
                    if(file != null) {
                        if(file.endsWith(".ms2") || file.endsWith(".mzXML")) {
                            String ms2filename = f + File.separator + file;
                            ms2files.add(ms2filename);
                            //System.out.println("adding " + ms2filename);
                        }
                    }
                }
            } else {
                ms2files.add(f);
            }        
        }
        return ms2files;
    }

    public static ArrayList<Double> getAllPrecursorMasses(ArrayList<String> ms2files, SearchParams params) throws Exception {
        ArrayList<Double> allprcmass = new ArrayList(20000);


        int numpeaksmin = params.getMinNumSpectra(); 
        int numpeaksmax = params.getMaxNumSpectra(); 
        double minprccharge = params.getMinPrecursorCharge();
        double maxprccharge = params.getMaxPrecursorCharge();
        double minprcmass = params.getMinPrecursorMass(); 
        double maxprcmass = params.getMaxPrecursorMass();

        for(Iterator<String> it = ms2files.iterator(); it.hasNext();) {
            String msFileName = it.next();
            ArrayList<PeakList> spectra = getSpectra(msFileName);
            Iterator<PeakList> itpl = spectra.iterator();
            
            while(itpl.hasNext()) {
                PeakList peaks = itpl.next();
                int numPeaks = peaks.numPeaks();

                if(numPeaks >= numpeaksmin && numPeaks <= numpeaksmax) {

                    for(Iterator<Zline> zlines = peaks.getZlines(); zlines.hasNext();) {
                        Zline z = zlines.next();
                        double precMass = z.getM2z();
                        int chargeState = z.getChargeState(); 
                        //int numPeaks = peaks.numPeaks();
                        if(chargeState < minprccharge || chargeState > maxprccharge) {
                            continue;
                        }

                        if(precMass >= minprcmass && precMass <= maxprcmass) {
                            allprcmass.add(new Double(precMass));                    
                        }
                    }
                }
            }
        }
        return allprcmass;
    }

    public static void main(String args[]) throws Exception {
        try {
            int numcpu = ManagementFactory.getOperatingSystemMXBean().getAvailableProcessors();
            //int numthread = numcpu;
            //int numthread = 4;
            int numthread = 1;

            String msFileName = args[0];
            String searchxml = "search.xml";
            if(args.length > 1) {
                searchxml = args[1];
            }
            if(args.length > 2) {
                numthread = Integer.parseInt(args[2]);
            }
            if(args.length >3)
            {
                PeptideHit.B_ION_FACTOR = Double.parseDouble(args[3]);
            }
            numthread = numthread > numcpu? numcpu : numthread;

            logbuffer.append("Number of thread\t" + numthread + "\n");
            SearchParams sp = new SearchParams(searchxml);
            SearchResult.setNumOutput(sp.getNumOutput());           

            ArrayList<String> ms2files = getMs2Files(msFileName);

            ArrayList<Double> allprecursormasses = getAllPrecursorMasses(ms2files, sp);
 
            ProteinDatabase pd = null;
            File pdindexfile = new File(sp.getDbName()+ ".pin");
            if(sp.isDatabaseIndexed() && pdindexfile.isFile()) {
                //pd = new IndexedProteinDatabase(sp.getDbName());
                System.out.println("Using indexed protein database");
                pd = new MemorizedIndexedProteinDatabase(sp.getDbName());
            }
            else if (sp.isUseProteoformDB())
            {
                System.out.println("Using proteoform database");
                pd = new ProteoformProteinDatabase(sp);
            }
            else {

                pd = new ProteinDatabase(sp.getDbName());
            }

            pd.setMaxMisCleavage(sp.getMaxInternalMisCleavage());
            pd.setProtease(sp.getProtease());
            pd.setEnzymeSpecificity(sp.getEnzymeSpecificity());

            long totalmemory = Runtime.getRuntime().totalMemory();
            long usedmemory  = totalmemory - Runtime.getRuntime().freeMemory();
            
            logbuffer.append("Total memory\t" + totalmemory + "\n");
            logbuffer.append("Memory used after loading the database\t" + usedmemory + "\n");
           
           
            // need to loop through all ms2 files 
            for(Iterator<String> msfileit = ms2files.iterator(); msfileit.hasNext();) {
                //ProteinDatabase pd = new IndexedProteinDatabase(args[0]);
                //ProlucidSearchEngine se = new ProlucidSearchEngine(sp, pd, msFileName);
                String ms2filename = msfileit.next();
                ProlucidSearchEngine se = new ProlucidSearchEngine(sp, pd, ms2filename);

                logbuffer.append("\nNow processing " + ms2filename + "\n");
                logbuffer.append("Total Memory after loading the search engine\t" + (Runtime.getRuntime().totalMemory() + "\tFree memory\t" + Runtime.getRuntime().freeMemory()) + "\n");
                logbuffer.append("Memory used after loading the search engine\t" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) + "\n");
                logbuffer.append("Host name\t" + getHostName() + "\n");
                long beginTime = System.nanoTime();
                TimeUtils timer = new TimeUtils();
                timer.startTiming();
                int numSpectraInFile = se.numSpectraInFile;
                logbuffer.append("TotalNumSpectra\t" + numSpectraInFile + "\n");

                logbuffer.append("NumSearched\t0\t" + numSpectraInFile + "\n" );
                System.out.print(logbuffer.toString()); 
                logbuffer = new StringBuffer(2000);

                ArrayList<Thread> threads = new ArrayList<Thread>();
                for(int i = 0; i < numthread; i++) {
                
                    Thread th = new Thread(new ProlucidThread(se)); 
                    threads.add(th);
                    th.start();
                } 
                for(Iterator<Thread> it = threads.iterator(); it.hasNext();) {
                    it.next().join();
                }


                long endTime = System.nanoTime();
                long timeUsed = (endTime - beginTime)/1000000;
                se.closeOutputFile();
                logbuffer.append("Total Memory after the search\t" + (Runtime.getRuntime().totalMemory() + "\tFree memory\t" + Runtime.getRuntime().freeMemory()) + "\n");
                logbuffer.append("Memory used after the search\t" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) + "\n");
                logbuffer.append("NumSpectraSearched: " + se.numSpectraSearched + "\t numSpectraIgnored: " + se.numSpectraIgnored + "\ttotalNumSpectra: " + se.numSpectraInFile+ "\n");
                logbuffer.append("NumResults " +se.nResults+"\t ignored results "+se.ignoredResults + "\n");
                logbuffer.append("Skipped because of precMass< minPrcMax: " +se.nMinPrcMass+"\t precMass>maxPrcMax: "+se.nMaxPrcMass +"\t pTrue>0.5: "+se.nPtrue +"\n");

                logbuffer.append("Time used for searching " + ms2filename + ": " + timeUsed + "\n");

                System.out.println(logbuffer.toString());
            }
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
            System.out.println("Problematic run on " + getHostName());
        }
    }
}


