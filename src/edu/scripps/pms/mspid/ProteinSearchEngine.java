/**
 * @file ProteinSearchEngine.java
 * This is the source file for edu.scripps.pms.mspid.ProteinSearchEngine
 * @author Tao Xu
 * @date $Date: 2009/11/04 00:45:24 $
 */
package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

import java.net.*;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.mspid.db.IndexedProteinDatabase;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.MzxmlSpectrumReader;
import org.jdom.JDOMException;
import edu.scripps.pms.util.io.TopDownSpectrumReader;

class ProteinSearchEngine {

    //public static final int DEFAULTNUMPEPTIDEHIT = 1000;
    //public static final int DEFAULTPEPTIDELENGTH = 40;
    private SearchParams params;
    private ProteinDatabase pdb;
    private MassCalculator mc;
    //private ProcessedPeakList ppl = null;
    public static final int HYPERGEOMETRYSCORE = 1;
    private static final String USAGE = "java ProteinSearchEngine ms2fileOrmzXMLfile";
    private String hostName;
    private long startTime;
    private long endTime;
    private double maxPrcMass;
    private double minPrcMass;

    private static double [] scan2mplush = new double[100000];
    private static int [] scan2charge = new int[100000];

    private static StringBuffer HEADER = new StringBuffer(1200);
    private static DistributionCalculator dc;
    public static final String HMLINE = "H\tMLINE\tPrimary_Rank\tSecondary_Rank\tTheoretical_Mass\tDelta_CN\t" +
            "Primary_Score(XCorr)\tSecondary_Score(ZScore)\tProbability_Score\tP_Score\tE_Score\tNumber_Of_Peaks_Matched\t" +
            "Number_Of_Peaks\tExtended_Sequence\n";

    public ProteinSearchEngine(SearchParams sp, ProteinDatabase db) { 
        this.params = sp;
        minPrcMass = params.getMinPrecursorMass(); 
        maxPrcMass = params.getMaxPrecursorMass();
        pdb = db;
        mc = new MassCalculator(sp);
        hostName = getHostName();
        getHeader();
    }

    private void getHeader() {
        HEADER.append("H\tSQTGenerator\tTopdown_ProLuCID\n");
        HEADER.append("H\tSQTGeneratorVersion\t0.3\n");
        HEADER.append("H\tComment ProLuCID is developed in the Yates laboratory at The Scripps Research Institute, La Jolla, CA\n");
        //HEADER.append("H\tComment ProLuCID ref. MOL CELL PROTEOMICS 5 (10): S174-S174 671 Suppl. S OCT 2006\n");
        //HEADER.append("H\tComment ProLuCID ref. Xu T, Venable JD, Park SK, Cociorva D, Lu B, Liao L, Wohlschlegel J, Hewel J, Yates JR 3rd\n");
        //HEADER.append("H\tComment Paralellization Program using PBS is submit_prolucid\n");
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
            HEADER.append("H\tStaticMod\t" + m.getResidue()+ "+" + m.getMassShift() + "\n");
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
                HEADER.append("H\tInternalDiffMod\t" + m.toString() + "\n");
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
        HEADER.append(HMLINE);

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
        for(Iterator<Zline> zlines = peaks.getZlines(); zlines.hasNext();) {
            Zline z = zlines.next();
            double precMass = z.getM2z();
            int chargeState = z.getChargeState(); 
            //int numPeaks = peaks.numPeaks();
            startTime = System.currentTimeMillis();
      
            ProcessedPeakList ppl = new ProcessedPeakList(peaks, z, params, mc); 
            if(precMass >= minPrcMass && precMass <= maxPrcMass && ppl.getPTrue() < 0.5) {
                //SearchResult sr = ppl.search(pdb);//search(peaks, params, z);
                SearchResult sr = pdb.search(ppl);//search(peaks, params, z);
         
            //System.out.println("TimeUsed: " + timeUsedInSeconds);
                sr.calcScores();
                sr.setHostName(hostName);
            //System.out.println(sr.outputResults());
                //endTime = System.currentTimeMillis();
                //long timeUsedInSeconds = (endTime - startTime);
                sr.setSearchTime(System.currentTimeMillis() - startTime);
                results.add(sr);
            }
            //System.out.println("ChargeState: " + chargeState + "\tPrecursorM2z: "
            //       + precMass + "\tNumber of peptideHits: " + sr.getNumPeptideHit());
       	}
        if(params.getChargeDisambiguation()) {
            results = chargeDisambiguation(results);
        }
        
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
    private static Iterator<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception  {
        Iterator<PeakList> peaklists = null;
        if(file.endsWith(".mzXML")) {
            MzxmlSpectrumReader sr = new MzxmlSpectrumReader(file);
            ArrayList<PeakList> spectra = new ArrayList<PeakList>(20000);
            for(Iterator<MzxmlPeakList> it = sr.getSpectra(2); it.hasNext();) {
                spectra.add(it.next()); 
            }
            peaklists = spectra.iterator();
            sr.closeDataFile();
            sr = null;
        } else {
            SpectrumReader sr = new SpectrumReader(file, "ms2");
            peaklists = sr.getSpectraList().iterator();
            sr.closeDataFile();
            sr = null;
        }
    
        return peaklists;
    }

    private static void readMs2(String ms2file) throws IOException, JDOMException, Exception {
        
        for (Iterator<PeakList> it = getSpectra(ms2file); it.hasNext();) {
            PeakList pl = it.next(); 
            int scan = pl.getLoscan();
            if(pl.getNumZlines() == 1) {
                double mh = 0;
                int charge = 0;
                for(Iterator<Zline> itz = pl.getZlines(); itz.hasNext();) {
                    Zline z = itz.next();
                    if(z.getChargeState() > 3) {
                        mh = z.getM2z();
                        charge = z.getChargeState();
                        System.out.println(scan + "\t" + charge + "\t"  + mh);
                    } else {
                        //System.out.println("charge state <= 3 " + charge);
                    }
                  
                }
                scan2mplush[scan] = mh;
                scan2charge[scan] = charge;
            }

        }
    }
    public static void main(String args[]) throws Exception {
        try {
            String msFileName = args[0];
            //scan2mplush = readMs2(msFileName);
            readMs2(msFileName);
            SearchParams sp = new SearchParams("search.xml");
            SearchResult.setNumOutput(sp.getNumOutput());           
 
            ProteinDatabase pd = null;
            if(sp.isDatabaseIndexed()) {
                pd = new IndexedProteinDatabase(sp.getDbName());
            } else {

                pd = new ProteinDatabase(sp.getDbName());
            }
            dc = new DistributionCalculator();
            //ProteinDatabase pd = new IndexedProteinDatabase(args[0]);
            ProteinSearchEngine se = new ProteinSearchEngine(sp, pd);
            //SpectrumReader sr = new SpectrumReader(msFileName, "ms2");
            //String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".ms2")) + ".sqt";
            //String outFileName =  "proteinid.txt";
            String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".")) + ".sqt";
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
            int numSpectraSearched = 0;
            int numSpectraIgnored = 0;
            int totalNumSpectra = 0;


            TopDownSpectrumReader sr = new TopDownSpectrumReader("xt/", "xls");
 
            for (Iterator<PeakList> it = sr.getSpectraList().iterator(); it.hasNext();) {
                i++;
                PeakList pl = it.next();        

                totalNumSpectra++;
                int numPeaks = pl.numPeaks();
                double mplush = scan2mplush[pl.getLoscan()];
                int charge = scan2charge[pl.getLoscan()];
        
                if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax && mplush != 0 ) {
                //timer.stopTiming(); 
                //System.out.println("Time used to read spectrum number " + i + ": " + timer.getTimeUsedMillis());
                //double prcMass = pl.getPrecursorMass();
                    pl.addZline(new Zline(charge, mplush));
                    List<SearchResult> results = se.search(pl);
                    
                    if(results.size() > 0) {
                        numSpectraSearched++;
                    } else {
                        numSpectraIgnored++;
                    }
                    for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {
                        
                        SearchResult rest = resultIt.next();
                        if(rest != null) { 
                            outFile.print(rest.outputResults());
                        }
                    }
               
                    System.gc();
                } else {
                    numSpectraIgnored++;
                } 
            }




            long endTime = System.currentTimeMillis();
            long timeUsed = endTime - beginTime;
            outFile.close();
            System.out.println("NumSpectraSearched: " + numSpectraSearched + "\t numSpectraIgnored: " + numSpectraIgnored + "\ttotalNumSpectra: " + totalNumSpectra);
            System.out.println("Time used for search: " + timeUsed);
            System.out.println();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
}


