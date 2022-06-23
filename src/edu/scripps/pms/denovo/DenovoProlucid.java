/**
 * @file DenovoProlucid.java
 * This is the source file for edu.scripps.pms.mspid.DenovoProlucid
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.denovo;

import edu.scripps.pms.util.seq.Fasta;

import edu.scripps.pms.mspid.*;
import java.net.*;
import java.util.*;
import java.io.*;
import java.lang.management.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.mspid.db.IndexedProteinDatabase;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.MzxmlSpectrumReader;
import org.jdom.JDOMException;

class DenovoProlucid {

    //public static final int DEFAULTNUMPEPTIDEHIT = 1000;
    //public static final int DEFAULTPEPTIDELENGTH = 40;
    private SearchParams params;
    private ProteinDatabase pdb;
    private MassCalculator mc;
    //private ProcessedPeakList ppl = null;
    public static final int HYPERGEOMETRYSCORE = 1;
    private static final String USAGE = "Usage: java DenovoProlucid ms2fileOrmzXMLfile";
    private String hostName;
    //private long startTime;
    private long endTime;
    private double maxPrcMass;
    private double minPrcMass;

    private static StringBuffer HEADER = new StringBuffer(1200);
    private static DistributionCalculator dc;
    private Iterator<PeakList> itr = null;
    private int numSpectraSearched = 0;
    private int numSpectraIgnored = 0;
    private PrintWriter outFile = null;
    private int numPeaksMin = 0; 
    private int numPeaksMax = 0; 
    private int numSpectraInFile = 0;
    protected static StringBuffer logbuffer = new StringBuffer(2000);
    
 
    public DenovoProlucid(SearchParams sp, ProteinDatabase db, String msFileName) 
                                                   throws IOException, JDOMException, Exception { 
        this.params = sp;
        minPrcMass = params.getMinPrecursorMass(); 
        maxPrcMass = params.getMaxPrecursorMass();
        pdb = db;
        mc = new MassCalculator(sp);
        hostName = getHostName();
        String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".")) + ".sqt";
        outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName)));
        getHeader();
        outFile.print(HEADER);
        ArrayList<PeakList> spectra = getSpectra(msFileName);
        numSpectraInFile =  spectra.size();
        itr = spectra.iterator();
        numPeaksMin = sp.getMinNumSpectra(); 
        numPeaksMax = sp.getMaxNumSpectra(); 
        if(dc == null) {
            dc = new DistributionCalculator();
        }

    }

    public PeakList getSpectrum() {
        if(itr != null && itr.hasNext()) {
            return itr.next(); 
        }
        return null; 
    }
    private void getHeader() {
        HEADER.append("H\tSQTGenerator\tProLuCID\n");
        HEADER.append("H\tSQTGeneratorVersion\t1.2.0\n");
        HEADER.append("H\tComment ProLuCID is developed in the Yates laboratory at The Scripps Research Institute, La Jolla, CA\n");
        HEADER.append("H\tComment ProLuCID ref. MOL CELL PROTEOMICS vol. 5(10): S174-S174 671 Suppl. S OCT 2006\n");
        HEADER.append("H\tComment ProLuCID ref. Xu T, Venable JD, Park SK, Cociorva D, Lu B, Liao L, Wohlschlegel J, Hewel J, Yates JR 3rd\n");
        HEADER.append("H\tComment Paralellization Program using PBS is submit_prolucid\n");
        HEADER.append("H\tComment Please send bug report or comments to Tao Xu by email taoxu@scripps.edu\n");
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
    }
    public static String getHostName() {
        String hostName = null;
        try {
            hostName  = InetAddress.getLocalHost().getHostName();
            logbuffer.append("host name " + hostName + "\n");
            if(hostName.indexOf(".") != -1) {
                hostName = hostName.substring(0, hostName.indexOf("."));
            }
        } catch (UnknownHostException e) {
            e.printStackTrace();
        }
        return hostName;

    }
    public void search(PeakList peaks) throws IOException {
        ArrayList<SearchResult> results = new ArrayList<SearchResult>();
        int numPeaks = peaks.numPeaks();
        if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax) {
            numSpectraSearched++;
        
        } else {
            numSpectraIgnored++;
            return;
        } 


        //System.out.print("NumSearched\t" + (numSpectraSearched+numSpectraIgnored) + "\t" + numSpectraInFile + "\t" );
        //PeptideHit.setFragMasses(mc.getFragMasses()); makes no big difference
        int totalNumProcessedSpectra = numSpectraSearched+numSpectraIgnored; 

        String logline = "NumSearched\t" + totalNumProcessedSpectra + "\t" + numSpectraInFile + "\t" ;
        logbuffer.append(logline);
        if(totalNumProcessedSpectra%20 == 0) {

            //System.out.print("NumSearched\t" + totalNumSpectra + "\t" + numSpectraInFile + "\t" );
            System.out.print(logbuffer.toString()); 
            logbuffer = new StringBuffer(2000);
        } 
        for(Iterator<Zline> zlines = peaks.getZlines(); zlines.hasNext();) {
            Zline z = zlines.next();
            double precMass = z.getM2z();
            int chargeState = z.getChargeState(); 
            //int numPeaks = peaks.numPeaks();
            if(chargeState < params.getMinPrecursorCharge() || chargeState > params.getMaxPrecursorCharge()) {
                continue;
            }
      
            ProcessedPeakList ppl = new ProcessedPeakList(peaks, z, params, mc); 
            logbuffer.append("FNumPeaks: " + ppl.getFinalNumPeaks() +
                             "\tScan#: " + peaks.getLoscan() + "\t+" +
                              chargeState + "\tNumPeaks: " + peaks.numPeaks() + 
                              "\tm2z: " + precMass + "\tptrue: " + ppl.getPTrue() + "\n");

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
            }
            //System.out.println("ChargeState: " + chargeState + "\tPrecursorM2z: "
            //       + precMass + "\tNumber of peptideHits: " + sr.getNumPeptideHit());
       	}
        if(params.getChargeDisambiguation()) {
            results = chargeDisambiguation(results);
        }
        
        outputSearchResults(results);
        //System.out.println();
        //return results;
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
    private static ArrayList<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception  {
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
    private void outputSearchResults(List<SearchResult> results) {
        
        synchronized(outFile) {
            for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {
                
                SearchResult rest = resultIt.next();
                if(rest != null) { 
                    outFile.print(rest.outputResults());
                }
                //outFile.flush();
            }
        }

    }
    public static void main(String args[]) throws Exception {
        try {
            int numcpu = ManagementFactory.getOperatingSystemMXBean().getAvailableProcessors();
            int numthread = numcpu;

            String msFileName = args[0];
            String searchxml = "search.xml";
            if(args.length > 1) {
                searchxml = args[1];
            }
            SearchParams sp = new SearchParams(searchxml);
            SearchResult.setNumOutput(sp.getNumOutput());           
 
            ProteinDatabase pd = null;
            if(sp.isDatabaseIndexed()) {
                pd = new IndexedProteinDatabase(sp.getDbName());
            } else {

                pd = new ProteinDatabase(sp.getDbName());
            }
            //ProteinDatabase pd = new IndexedProteinDatabase(args[0]);
            DenovoProlucid se = new DenovoProlucid(sp, pd, msFileName);

            logbuffer.append("Number of thread\t" + numthread + "\n");
            long beginTime = System.nanoTime();
            TimeUtils timer = new TimeUtils();
            timer.startTiming();
            logbuffer.append("TotalNumSpectra\t" + se.numSpectraInFile + "\n");

            ArrayList<Thread> threads = new ArrayList<Thread>();
            for(int i = 0; i < numthread; i++) {
               
                Thread th = new Thread(new DenovoProlucidThread(se)); 
                th.start();
                threads.add(th);
            } 
            for(Iterator<Thread> it = threads.iterator(); it.hasNext();) {
                it.next().join();
            }


            long endTime = System.nanoTime();
            long timeUsed = (endTime - beginTime)/1000000;
            se.closeOutputFile();
            logbuffer.append("NumSpectraSearched: " + se.numSpectraSearched + "\t numSpectraIgnored: " + se.numSpectraIgnored + "\ttotalNumSpectra: " + se.numSpectraInFile + "\n");
            logbuffer.append("Time used for search: " + timeUsed + "\n");
            System.out.println(logbuffer.toString());
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
            System.out.println("Problematic run on " + getHostName());
        }
    }
}


