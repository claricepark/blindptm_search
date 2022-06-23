package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.proteoform.Proteoform;
import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSequenceWithPTMs;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.proteoform_dbindex.model.PTMCodeObj;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import gnu.trove.map.hash.TIntByteHashMap;
import gnu.trove.map.hash.TIntCharHashMap;
import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.apache.poi.ss.formula.functions.T;
import org.hsqldb.Index;
import org.jdom.JDOMException;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by Titus Jung titusj@scripps.edu on 9/12/19.
 */
public class ProteoformProteinDatabase extends ProteinDatabase {

    private ProteoformDBIndexInterface proteoformDBIndex;
    private ExtendedAssignMass extendedAssignMass;

    public ProteoformProteinDatabase(SearchParams sp) throws IOException {
        super(sp.getDatabaseName());

       /* List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
        loggers.add(LogManager.getRootLogger());
        for ( Logger logger : loggers ) {
            logger.setLevel(Level.OFF);
        }*/
       //Logger.getLogger("edu.scripps.yates.dbindex.DBIndexImpl").setLevel(Level.OFF);
        String directory = sp.getDatabaseName().substring(0,sp.getDatabaseName().lastIndexOf(File.separator));
        int maxMisclevages = sp.getMaxInternalMisCleavage();
        double minPrcMass = sp.getMinPrecursorMass();
        double maxPrcMass = sp.getMaxPrecursorMass();
        String residues = sp.getProtease().getResidues();
        int maxNumVariationsPerPeptide = sp.getMaxAlter();
        int offset = sp.getEnzyme().getType() == true? 1: 0;
        String uniprotVersion = sp.getUniprotVersion();
        DBIndexSearchParams params =  DBIndexImpl.getDefaultDBIndexParamsForProteoformAnalysis(sp.getDbName(),
                maxMisclevages, minPrcMass,maxPrcMass, residues, null,offset,false,
                uniprotVersion, null, directory);
        proteoformDBIndex = new ProteoformDBIndexInterface(params,maxNumVariationsPerPeptide);
        extendedAssignMass = ExtendedAssignMass.getInstance(false,null);
    }

    public SearchResult search(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        double [] precMasses = ppl.getMassCalculator().getPrecMasses();
        //isDeCharged = ppl.isDeCharged();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();
        //MassCalculator mc = ppl.getMassCalculator();
        //System.out.println("DEBUG><<>><>" +ppl.peakList.getLoscan());
        //double mass = prcMass;

        if(numIsotopes == 0) { // for low resolution, traditional sequest like
        //    double highLimit = (prcMass + params.getHighPrecursorTolerance()/1000);
        //    double lowLimit = (prcMass - params.getLowPrecursorTolerance()/1000);
            getPeptideHitsTol(sr, prcMass, params.getHighPrecursorTolerance()/1000);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data
            double diffs = prcMass*acc;
          //  double highLimit = prcMass + diffs;
          //  double lowLimit = prcMass - diffs;
            getPeptideHitsTol(sr, prcMass, diffs);
        } else { //for non deisotoped high resolution data
        //    double [] highLimits = new double[numIsotopes];
         //   double [] lowLimits = new double[numIsotopes];
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
               // lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13;
              //  highLimits[i] = lowLimits[i] + diffs;
                getPeptideHitsTol(sr, prcMass-i*MassSpecConstants.MASSDIFFC12C13, diffs );

            }

        }

        return sr;
    }

    public synchronized List<IndexedSequence> queryPeptides(double prcMass, double tolerance) throws DBIndexStoreException {
        return  proteoformDBIndex.getSequences(prcMass, tolerance);
    }

    public  void  getPeptideHitsTol(SearchResult sr, double prcMass, double tolerance)  {
        try {

            List<IndexedSequence> indexedSequenceList = queryPeptides(prcMass,tolerance);
            for(IndexedSequence iseq: indexedSequenceList)
            {
               int fastaID = iseq.getProteinIds().get(0);
               IndexedProtein protein=  proteoformDBIndex.getIndexedProteinById(fastaID);
               String left = "";
               String right = "";
               int start =0;
               if(!iseq.getResLeft().equals("---"))
               {
                   start =3;
                   left = iseq.getResLeft();
               }
               if(!iseq.getResRight().equals("---"))
               {
                   right = iseq.getResRight();
               }
               Fasta fasta = new Fasta(protein.getFastaDefLine(),  left+ iseq.getSequence() +right);
            //    System.out.println("<<<>>> "+iseq.getSequence());
               int end = start + iseq.getSequenceLen()-1;
                if(iseq.isIsModified())
                {
                 //  System.out.println("adding >><<> "+iseq.getModSequence());
                    IndexedSequenceWithPTMs iseqPTM = (IndexedSequenceWithPTMs) iseq;
                    List<PTM> ptmSet = iseqPTM.getPtms();
                    Map<Integer,DiffMod> diffModMap = new HashMap<>();
                    Map<Integer,List<String>> polyMap = new HashMap<>();
                    for(PTM ptm : ptmSet)
                    {
                        int index  = ptm.getPosInPeptide()-1;
                        PTMCodeObj codeObj = extendedAssignMass.getPTMbyPTMCode(ptm.getPtmCode());
                      //  System.out.println("<<>> "+codeObj.getDescription());
                        if(codeObj.getPtmMassDiff()!=null)
                        {
                            if(index>=iseq.getSequence().length())
                            {
                                System.out.println(iseq.getModSequence() + "\t" + protein.getFastaDefLine());

                                System.out.println("DEBUG>><<<>"+sr.getProcessedPeakList().getPeakList().getLoscan());
                            }
                            else
                            {
                                char c = iseq.getSequence().charAt(index);
                                double mass = codeObj.getPtmMassDiff();
                                //System.out.println("<<<c "+c+" "+mass+" "+ptm.getPosInPeptide());
                                DiffMod diffMod = new DiffMod(mass, c);
                                diffModMap.put(index,diffMod);
                            }

                        }
                        else
                        {
                         //   System.out.println("<<<c "+codeObj.getDescription()+" "+ptm.getPosInPeptide());

                            String poly = codeObj.getDescription();
                           // System.out.println(">> "+poly);
                          //  System.out.println(">>>>"+iseq.getSequence());
                          //  System.out.println(">>> "+iseq.getModSequence());
                            List<String> polyList = polyMap.getOrDefault(index,new ArrayList<String>());
                            polyList.add(poly);
                            polyMap.put(index,polyList);
                        }
                    }
                  //  List<Integer> keyList = new ArrayList<>(polyMap.keySet());
                   // List<List<Integer>> subSetList = getAllSubsets(keyList);`
                   /* for(List<Integer> subset : subSetList)
                    {*/
                        ModifiedPeptideHit mph = new ModifiedPeptideHit(fasta,start,start+iseq.getSequenceLen()-1,new Modifications());

                        for(Map.Entry<Integer,DiffMod> entry: diffModMap.entrySet())
                        {
                            mph.setDiffMod(entry.getKey(), entry.getValue());
                        }
                        mph.setPolyMorphMap(polyMap);
                        mph.setIseq(iseqPTM);
                        sr.addMPeptideHit(mph);
                   // }
                }
                else
                {
                    sr.addPeptideHit(fasta, start, end);
                }
            }
        } catch (DBIndexStoreException e) {
            e.printStackTrace();
        }
    }

    public static List<List<Integer>> getAllSubsets(List<Integer> input)
    {
        List<List<Integer>> result = new ArrayList<>();
        int n = result.size();
        for (int i = 0; i < (1<<n); i++)
        {
            List<Integer> subset = new ArrayList<>();
            for (int j = 0; j < n; j++)
            {
                if ((i & (1 << j)) > 0)
                {

                    subset.add(input.get(j));
                }
            }
            result.add(subset);
        }

        return result;
    }

    public static void main(String [] args) throws IOException, JDOMException {
        String searchXMLPath = args[0];
        SearchParams sp = new SearchParams(searchXMLPath);
        ProteoformProteinDatabase db = new ProteoformProteinDatabase(sp);
    }


}
