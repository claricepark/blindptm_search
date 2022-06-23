package edu.scripps.pms.protinf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.Protein;
//import edu.scripps.pms.util.spectrum.PeakList;
//import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.seq.PrefixDb;
import java.util.logging.Level;
import java.util.logging.Logger;
//import edu.scripps.pms.mspid.MassCalculator;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class ProteinInference {


    private static String outputPath = ".";

    public static String HELP = "proteininferencer result_folder_list_file protein_database_file [peptdie_FDR]";
    // C_elegan
    //public static final String database = "/lustre/people/applications/yates/dbase/wormbase_c-elegans_225_NoGeneDuplicates_05-10-2011_reversed.fasta";
    public static String database = "/lustre/people/applications/yates/dbase/UniProt_mouse_05-17-2011_reversed.fasta";
    //public static final String database = "/lustre/people/applications/yates/dbase/UniProt_Human_05-17-2011_reversed.fasta";
    //public static final String database = "/data/2/rpark/ip2_data/taoxu/database///17PM_pombe030305.fasta";
    //public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/protinf/celegan_mitochodira.txt";
    public static String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/protinf/MouseMitochodrialV1.txt";
    //public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/protinf/test.txt";
//    public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/protinf/protinf_test.txt";
    public static int peptideScoreType = 1;
    //public static HashMap<String, Fasta> ac2Fasta = new HashMap();
    public static ArrayList<String> folders = new ArrayList(); 
    //private static MassCalculator mc = new MassCalculator();
    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println("usage: java ProteinInference file_contains_all_the_folders output_file_name");
	System.exit(-1);
    }

    public static int getNumExperimentsAdded() {
        return folders.size();
    }

    public static String getHeaderByExperiment(String hd) {
        StringBuffer sb = new StringBuffer();
        for(int i = 0 ; i < folders.size(); i++) {
            sb.append(hd);
            sb.append(folders.get(i));
            sb.append("\t");
        }
        return sb.toString();
    }

    public static void assignPeptideItemDecoyStatus(ArrayList<ProteinGroup> pgs) {
        for(Iterator<ProteinGroup> it = pgs.iterator(); it.hasNext();) {
            ProteinGroup pg = it.next();
           
            for(Iterator<PeptideItem> pepit = pg.getPeptideItems().iterator(); pepit.hasNext();) {
                PeptideItem pi = pepit.next();
                pi.isReverseHit(pg.isReverseHit());
            }
        }
    }
    public static int getPeptideLocation(String peptide, Fasta protein) {
        return protein.getSequence().indexOf(peptide);
    }

    private static String getCummulativeProteinCounts(ArrayList<ProteinGroup> pgs) {
        StringBuffer sb = new StringBuffer();
        StringBuffer sb0 = new StringBuffer(); sb0.append("Protein FDR 0\t");
        StringBuffer sb001 = new StringBuffer(); sb001.append("Protein FDR 0.001\t");
        StringBuffer sb01 = new StringBuffer(); sb01.append("Protein FDR 0.01\t");
        StringBuffer sb02 = new StringBuffer(); sb02.append("Protein FDR 0.02\t");
        StringBuffer sb03 = new StringBuffer(); sb03.append("Protein FDR 0.03\t");
        StringBuffer sb05 = new StringBuffer(); sb05.append("Protein FDR 0.05\t");
        StringBuffer sb10 = new StringBuffer(); sb10.append("Protein FDR 0.1\t");
  
        HashSet<ProteinGroup> conf0 = new HashSet();
        HashSet<ProteinGroup> conf001 = new HashSet();
        HashSet<ProteinGroup> conf01 = new HashSet();
        HashSet<ProteinGroup> conf02 = new HashSet();
        HashSet<ProteinGroup> conf03 = new HashSet();
        HashSet<ProteinGroup> conf05 = new HashSet();
        HashSet<ProteinGroup> conf10 = new HashSet();
        for(int i = 0; i < folders.size(); i++) {
            for(Iterator<ProteinGroup> it = pgs.iterator(); it.hasNext();) {
               ProteinGroup pg = it.next();
               double fdr = pg.setFinalGlobalFDr();
               
               if(pg.getRepresentative().getPeptidesByExperiment()[i].size() > 0 && !pg.isReverseHit()) { // identified in this experiment
       
                    if(fdr <= 0) {
                        conf0.add(pg);
                    } else if(fdr <= 0.001) {
                        conf001.add(pg);
                    } else if(fdr <= 0.01) {
                        conf01.add(pg);
                    } else if(fdr <= 0.02) {
                        conf02.add(pg);
                    
                    } else if(fdr <= 0.03) {
                        conf03.add(pg);
                    
                    } else if(fdr <= 0.05) {
                        conf05.add(pg);
                    
                    } else if(fdr <= 0.1) {
                        conf10.add(pg);
                    } 
               } 
            }
            HashSet<ProteinGroup> cpgs = new HashSet();
            cpgs.addAll(conf0);
            sb0.append(cpgs.size() + "\t");
            cpgs.addAll(conf001);
            sb001.append(cpgs.size() + "\t");
            cpgs.addAll(conf01);
            sb01.append(cpgs.size() + "\t");
            cpgs.addAll(conf02);
            sb02.append(cpgs.size() + "\t");
            cpgs.addAll(conf03);
            sb03.append(cpgs.size() + "\t");
            cpgs.addAll(conf05);
            sb05.append(cpgs.size() + "\t");
            cpgs.addAll(conf10);
            sb10.append(cpgs.size() + "\t");
            
         
        }
        sb.append("Protein FDR\t" +  getHeaderByExperiment("") + "\n");
        sb.append(sb0.toString() + "\n");  
        sb.append(sb001.toString() + "\n");  
        sb.append(sb01.toString() + "\n");  
        sb.append(sb02.toString() + "\n");  
        sb.append(sb03.toString() + "\n");  
        sb.append(sb05.toString() + "\n");  
        sb.append(sb10.toString() + "\n");  
        return sb.toString();
    }
    public static HashMap<String, Fasta> getProteinMap(String database) throws IOException {
        HashMap<String, Fasta> ac2fasta = new HashMap(1000000);
        FileInputStream fis = new FileInputStream(new File(database));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        while(fastas.hasNext()) {
            Fasta f = fastas.next();
            String myac = f.getAccession();
   //System.out.println("protien ac " + myac + "\tdescription " + f.getDescription());
            //if(myac != null && !myac.startsWith("Reverse")) {
            if(myac != null) {
                ac2fasta.put(myac, f);
            }
        }
        return ac2fasta;
    }

    public static ArrayList<ProteinGroup> groupProteins(HashSet<ProteinItem>  pis) {
//System.out.println("Number of ProteinItems before grouping " + pis.size());
        ArrayList<ProteinGroup> pgs = new ArrayList(100000);
        ArrayList<ProteinItem> temppis = new ArrayList();
        temppis.addAll(pis);
        Collections.sort(temppis);
        int numpgs = 0;
        int numreversepgs = 0;
//System.out.println("GroupProtein Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tsumZScore\tIsReverseHit");
        for(int i = 0; i < temppis.size(); i++){
            ProteinItem current = temppis.get(i);
            if(current != null) {
                ProteinGroup pg = new ProteinGroup();
                numpgs++;
                pg.addProteinItem(current);
                if(pg.isReverseHit()) {
                    numreversepgs++;
                }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
                pgs.add(pg);
                temppis.set(i, null);
                HashSet<PeptideItem> currentset = current.getPeptideItems();
                for(int j = i+1; j < temppis.size(); j++) {
                
                    ProteinItem next = temppis.get(j);
                    if(next != null) {

                        HashSet<PeptideItem> nextset = next.getPeptideItems();
                        if(currentset.size() == nextset.size() && currentset.containsAll(nextset)) {
                            pg.addProteinItem(next);
                            temppis.set(j, null);
                        } 

                    }
                }
            }
        }        

System.out.println("Final number of protein groups added " + numpgs);
        return pgs;
    }

    private static void outputProteinConfidence(ArrayList<ProteinGroup> allpgs) { 
        HashMap<PeptideItem, HashSet<ProteinGroup>> peptide2pgs = new HashMap();
        
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            ProteinItem pi = pg.getRepresentative();
            Iterator<PeptideItem> it = pi.getPeptideItems().iterator();
            while(it.hasNext()) {
                PeptideItem pepitem = it.next();
                HashSet<ProteinGroup> pgs = peptide2pgs.get(pepitem);
                if(pgs == null) {
                    pgs = new HashSet();
                }
                pgs.add(pg);
            }
            
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
        }
        Iterator<PeptideItem> it = peptide2pgs.keySet().iterator();
        while(it.hasNext()) {
            PeptideItem pi = it.next();
            HashSet<ProteinGroup> pgs = peptide2pgs.get(pi);
            
        } 
    }

    // the allpgs should be sorted when this function is called, scoretype 2 for sumZScore, 5 for averageZScore, 6 for confidenceSum 
    // 7 for confidenceProduct 
    private static void outputProteinConfidence(ArrayList<ProteinGroup> allpgs, int scoretype) throws IOException {
       
        double [] localfdr = new double[allpgs.size()];
        int [] globlefdr = new int[allpgs.size()];
        boolean [] reversehits = new boolean[allpgs.size()];
        int numreversepgs = 0;
        int numpgs = 0;
        PrintWriter protFile = new PrintWriter(new BufferedWriter(new FileWriter(outputPath + File.separator+"protinf_proteingroups.txt")));
//System.out.println("Number of Protein Groups Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tIsReverseHit");
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            numpgs++;
            if(pg.isReverseHit()) {
                reversehits[i] = true;
                globlefdr[i] = ++numreversepgs;
                
            }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
        }
        int numcloseby = 100; 
        for(int i = 0; i < reversehits.length; i++) {
            int temprev = 0;
            for(int j = i - numcloseby; j < i+numcloseby; j++) {
                if(j > -1 && j < reversehits.length && reversehits[j]) {
                    temprev++;
                }
            }
            localfdr[i] = temprev;
        }

       for(int i = 0; i < localfdr.length; i++) {
           int before = i >= numcloseby? numcloseby : i;
           int after = localfdr.length - i >= numcloseby? numcloseby : (localfdr.length - i);
//System.out.print(i + "\t" + localfdr[i] + "\t");
           localfdr[i] = (localfdr[i]+0.0)/(1 + before + after);
           if(i > 0) {
//System.out.println(localfdr[i] + "\t" + localfdr[i-1]); 
               if(localfdr[i] < localfdr[i-1]) {       
                   localfdr[i] = localfdr[i-1];
               }
           }
       }
      
//System.out.println("Protein Group ID\tNumber of Members\tOccurrence\tSorting Code\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tNum Identified Tryptic Peptides\tNum Typtic Peptides\tProtein Length\tIsReverseHit\tRepresentative Protein accession\tDefline\tAll accessions\tAll Deflines");
//protFile.println("Protein Group ID\tNumber of Members\tOccurrence\tSorting Code\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tNum Identified Tryptic Peptides\tNum Typtic Peptides\tProtein Length\tIsReverseHit\tRepresentative Protein Accession\tRepresentative Protein Gene Name\tDefline\tAll accessions\tAll Deflines\t" + getSpectrumCountHeaderByExperiment() + getPeptideCountHeaderByExperiment());
protFile.println("Protein Group ID\tNumber of Members\tOccurrence\tSorting Code\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tNum Unique Peptide Sequences\tNum Identified Tryptic Peptides\tNum Typtic Peptides\tProtein Length\tRepresentative Protein Sequence Coverage (%)\tIsReverseHit\tRepresentative Protein Accession\tRepresentative Protein Gene Name\tDefline\tAll accessions\tAll Deflines\t" + getHeaderByExperiment("Spectrum_Count_") + getHeaderByExperiment("Peptide_Count_"));
        StringBuffer summarysb = new StringBuffer();
        StringBuffer proteinlistsb = new StringBuffer();
       
        int numgfr0 = 0;
        int numgfr01 = 0;        
        int numgfr1 = 0;        
        int numgfr2 = 0;        
        int numgfr3 = 0;        
        int numgfr5 = 0;        
        int numgfr10 = 0;        
 
        int numforward = 0;
        int numreverse = 0; 
        int pgid = 0;  // for protein group id
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            if(pg.isReverseHit()) {
                numreverse++; 
            } else {

                numforward++;
            }
    // the allpgs should be sorted when this function is called, scoretype 2 for sumZScore, 5 for averageZScore, 6 for confidenceSum 
    // 7 for confidenceProduct 
            switch(scoretype) {
                case 5 : pg.setAverageZScoreConfidence(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 2 : pg.setSumZScoreConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
                default : pg.setSumZScoreConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
            }

            float seqcov =  pg.getRepresentative().getSeqCov();
            String repacc = pg.getRepresentativeAcc();
            String repdefline = pg.getRepresentativeDefline();
            String accs = pg.getAccs();
            String deflines = pg.getDeflines();
            int numidedtyppeps = pg.getRepresentative().getNumIdentifiedTrypticPeptides(); 
            int numyppeps = pg.getRepresentative().getNumTrypticPeptides(); 
            int numuniquepepseq = pg.getRepresentative().getNumUniquePeptideSequence(); 
            double gfdr = (numreverse/(0.0 + numforward));
            double lfdr = (1-localfdr[i]);
            String repgenename = pg.getRepresentative().getFasta().getGeneName();
//System.out.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + numforward + "\t" + numreverse + "\t" + (numreverse/(0.0 + numforward)) + "\t" + (1-localfdr[i]) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore()
//        + "\t" + numidedtyppeps + "\t" + numyppeps + "\t" + pg.getRepresentativeLength() + "\t" + pg.isReverseHit() + "\t" + repacc + "\t" + repdefline + "\t" + accs + "\t" + deflines);
protFile.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + pg.getOccurrenceAndSortingCode() + "\t" + numforward + "\t" + numreverse + "\t" + gfdr + "\t" + lfdr + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore()
        + "\t" + numuniquepepseq + "\t" + numidedtyppeps + "\t" + numyppeps + "\t" + pg.getRepresentativeLength() + "\t" + seqcov + "\t" + pg.isReverseHit() +  "\t" + repacc + "\t" + repgenename + "\t" + repdefline + "\t" + accs + "\t" + deflines + 
        "\t" + pg.getRepresentative().getSpectrumCountInfo() + pg.getRepresentative().getPeptideCountInfo());

            pg.setFinalGlobalFDr(gfdr);
            pg.setFinalConfidence(lfdr);
            pg.setProteinGroupId(i+1);

            numgfr0 = gfdr <= 0 ? numforward : numgfr0;  
            numgfr01 = gfdr <= 0.001 ? numforward : numgfr01;  
            numgfr1 = gfdr <= 0.01 ? numforward : numgfr1;  
            numgfr2 = gfdr <= 0.02 ? numforward : numgfr2;  
            numgfr3 = gfdr <= 0.03 ? numforward : numgfr3;  
            numgfr5 = gfdr <= 0.05 ? numforward : numgfr5;  
            numgfr10 = gfdr <= 0.1 ? numforward : numgfr10;  
        }
        System.out.println("\n\n\nGlobal FDR\t0\t0.001\t0.01\t0.02\t0.03\t0.05\t0.1");
        System.out.println("Number of Forward Protein Groups\t" + numgfr0 +  "\t" + numgfr01 + "\t" + numgfr1 + "\t" + numgfr2 + "\t" + numgfr3 + "\t" + numgfr5 + "\t" + numgfr10);
        protFile.close(); 

    }
    private static String getProteinInfo(ArrayList<ProteinGroup> pgs) {
        StringBuffer allsb = new StringBuffer();
        StringBuffer repsb = new StringBuffer();
        int numpgs = pgs.size();
        double protconf = 0;
        double protfdr = 0;
        for(int i = 0; i < pgs.size(); i++) {
            ProteinGroup pg = pgs.get(i);
            repsb.append(pg.getRepresentativeAcc());
            repsb.append("|||");
            allsb.append(pg.getAccs());
            allsb.append("|||"); 
            if(i == 0) {
                protfdr = pg.setFinalGlobalFDr();
                protconf = pg.getFinalConfidence();

            }
        }
        return protconf + "\t" + protfdr  + "\t" + numpgs + "\t" + repsb.toString() + "\t" + allsb.toString() + "\t";
    }
    private static double getProteinFdr(ArrayList<ProteinGroup> pgs) {
        int numpgs = pgs.size();
        for(int i = 0; i < pgs.size(); i++) {
            ProteinGroup pg = pgs.get(i);
            return pg.setFinalGlobalFDr();

        }
        return 0;
    }
    private static int getTotalSpectrumCount(HashSet<PeptideItem> pis) {
        int totalspectrumcount = 0;
        for(Iterator<PeptideItem> it = pis.iterator(); it.hasNext();) {
            totalspectrumcount += it.next().getTotalSpectrumCount();
        }
        return totalspectrumcount;
    }
    private static int getNumPeptideSequences(HashSet<PeptideItem> pis) {
//System.out.println("PeptideItem number: " + pis.size());
        HashSet<String> pepseqs = new HashSet();
        for(Iterator<PeptideItem> it = pis.iterator(); it.hasNext();) {
            PeptideItem pi = it.next();
            pepseqs.add(pi.getSequence()); // or getSeqWithoutMod() ?
            //pepseqs.add(pi.getSeqWithoutMod()); // or getSeqWithoutMod() ?
//System.out.println("pepseq: " + pi.getSequence());
        }
//System.out.println("PeptideSequence number: " + pepseqs.size());
        return pepseqs.size();
    }
    // [0] and [1] are the number of forward and decoy PeptideItems,  [2] and [3] are the number of forward and decoy peptide sequences and [4] and [5] are number of forward and decoy spectrum counds
    private static int [] countNumForwardAndDecoyHits (HashSet<PeptideItem> pis) {
        int [] forwardanddecoyhits = new int[6];
        HashSet<String> forwardseqs = new HashSet();
        HashSet<String> decoyseqs = new HashSet();
        for(Iterator<PeptideItem> it = pis.iterator(); it.hasNext(); ) {
            PeptideItem pi = it.next();
            if(pi.isReverseHit()) {
                forwardanddecoyhits[1]++;
                decoyseqs.add(pi.getSequence());
                forwardanddecoyhits[5] += pi.getTotalSpectrumCount();
                
            } else {
                forwardanddecoyhits[0]++;
                forwardseqs.add(pi.getSequence());
      
                forwardanddecoyhits[4] += pi.getTotalSpectrumCount();
            }
        } 
        forwardanddecoyhits[2] = forwardseqs.size(); 
        forwardanddecoyhits[3] = decoyseqs.size(); 
        return forwardanddecoyhits;
    } 
    public static String getSpectrumCountByExperiment(ArrayList<PeptideItem> pis) {
        HashSet<PeptideItem> pitems = new HashSet();
        pitems.addAll(pis);
        int [] spectrumcounts = new int[getNumExperimentsAdded()]; 
        int totalspc = 0;
        for(Iterator<PeptideItem> it = pitems.iterator(); it.hasNext();) {
            PeptideItem pi = it.next();
            Peptide [] pipts = pi.getPeptides();
            for(int i = 0; i < pipts.length; i++) {
                Peptide p = pipts[i];
                if(p != null) {
                    int spc =  p.getRedundancyValue();
                    spectrumcounts[i] += spc;
                    totalspc += spc;
                }
            }
    
        }        

        StringBuffer sb = new StringBuffer();
        sb.append(pitems.size() + "\t");
        sb.append(totalspc + "\t");
        for(int i = 0; i < spectrumcounts.length; i++) {
            sb.append(spectrumcounts[i] + "\t");
        }
        return sb.toString();
       
   
    }
    private static void outputPeptideWithProteinConfidence(ArrayList<ProteinGroup> allpgs) throws IOException {
        HashMap<PeptideItem, ArrayList<ProteinGroup>> pep2pgs = new HashMap();
        PrintWriter pepFile = new PrintWriter(new BufferedWriter(new FileWriter(outputPath + File.separator+ "protinf_peptide.txt")));
        //pepFile.println("Peptide Sequence\tBest Hit Charge State\tTryptic Status\tIs Decoy Hit\tModification Status\tSorting Code\tOccurrence\tTotal Spectrum Count\tPeptide Confidence\tProtein Confidence\tProtein FDR\tNum Protein Groups\tProtein Group Representatives\tProtein Members\t" + getHeaderByExperiment("Spectrum_Count_") );
        pepFile.println("Peptide Sequence\tBest Hit Charge State\tTryptic Status\tIs Decoy Hit\tModification Status\tSorting Code\tOccurrence\tBest Hit XCorr\tBest Hit DeltaCN\tBest Hit Sp\tBest Hit Delat Mass (ppm)\tBest Hit Scan\tTotal Best Hit Spectrum Count\tPeptide Confidence\tProtein Confidence\tProtein FDR\tNum Protein Groups\tProtein Group Representatives\tProtein Members\tNum Charge State\tTotal Spectrum Count\t" + getHeaderByExperiment("Spectrum_Count_") );
        HashSet<PeptideItem> allpeptideitems = new HashSet();
// return seq + "\t" + charge + "\t" + getSortingCode() + "\t" + getOccurrence() + "\t" + getTotalSpectrumCount() + "\t" + getConfidence();
        HashMap<String, PeptideItem> pepseq2bestpeptideitem = new HashMap(); // to ensure only output one peptide per peptide sequence
        HashMap<String, ArrayList<PeptideItem>> pepseq2peptideitems = new HashMap(); // to ensure output all spectrum counts for each peptide sequence
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            for(Iterator<PeptideItem> it = pg.getRepresentative().getPeptideItems().iterator(); it.hasNext();) {
               PeptideItem pi = it.next();
               allpeptideitems.add(pi);
               String pepseq = pi.getSequence();
               ArrayList<PeptideItem> pitems = pepseq2peptideitems.get(pepseq);
               if(pitems == null) {
                   pitems = new ArrayList();
                   pepseq2peptideitems.put(pepseq, pitems);
               } 
               pitems.add(pi);  // these are for spectrum counting, needs all PeptideItem, not just the best

                               
               // the following is for the best PeptideItem
               PeptideItem bestpi = pepseq2bestpeptideitem.get(pepseq);
               if(bestpi == null) {
                   bestpi = pi;
               } else {
                   if(pi.getConfidence() > bestpi.getConfidence()) {
                       bestpi = pi;
                   } else if(pi.getConfidence() == bestpi.getConfidence()) {
                       //prefer charge 2 over charge 3 and charge 1
                       int bestpicharge = bestpi.getCharge();
                       int picharge = pi.getCharge();
                       if(picharge == 2) {
                           bestpi = pi;
                       } else if(bestpicharge == 1) {
                           bestpi = pi;
                       }
                      
                   }
               } 
               pepseq2bestpeptideitem.put(pepseq, bestpi);
               ArrayList<ProteinGroup> pgs = pep2pgs.get(pi);
               if(pgs == null) {
                   pgs = new ArrayList();
                   pep2pgs.put(pi, pgs);
               } 
               pgs.add(pg);
            }
        }


        //Set<PeptideItem> pis = pep2pgs.keySet();
        Collection<PeptideItem> pis = pepseq2bestpeptideitem.values();
        for(Iterator<PeptideItem> it = pis.iterator(); it.hasNext();) {
            PeptideItem pi = it.next();
            String pepseq = pi.getSequence();
            String pepinfo = pi.getPeptideInfo();
            //String spectrumcountinfo = pi.getSpectrumCountByExperiment();
            String spectrumcountinfo = getSpectrumCountByExperiment(pepseq2peptideitems.get(pepseq));
            pepFile.print(pepinfo);
            ArrayList<ProteinGroup> pgs = pep2pgs.get(pi);
            pepFile.print(getProteinInfo(pgs));
            pepFile.println(spectrumcountinfo);
        
        }

        HashSet<PeptideItem> protConf0 = new HashSet();
        HashSet<PeptideItem> protConf001 = new HashSet();
        HashSet<PeptideItem> protConf01 = new HashSet();
        HashSet<PeptideItem> protConf02 = new HashSet();
        HashSet<PeptideItem> protConf03 = new HashSet();
        HashSet<PeptideItem> protConf05 = new HashSet();
        HashSet<PeptideItem> protConf10 = new HashSet();
        for(Iterator<PeptideItem> it = allpeptideitems.iterator(); it.hasNext();) {
            PeptideItem pi = it.next();
            String pepinfo = pi.getPeptideInfo();
            ArrayList<ProteinGroup> pgs = pep2pgs.get(pi);
        
            double protfdr = getProteinFdr(pgs);
            if(protfdr <= 0) {
                protConf0.add(pi);
            } else if(protfdr <= 0.001) {
                protConf001.add(pi);
            } else if(protfdr <= 0.01) {
                protConf01.add(pi);
            } else if(protfdr <= 0.02) {
                protConf02.add(pi);
            } else if(protfdr <= 0.03) {
                protConf03.add(pi);
            } else if(protfdr <= 0.05) {
                protConf05.add(pi);
            } else if(protfdr <= 0.1) {
                protConf10.add(pi);
            }
        }


        HashSet<PeptideItem> allpis = new HashSet();
        allpis.addAll(protConf0); 
        int [] numprotconf0 = countNumForwardAndDecoyHits(protConf0); 
        allpis.addAll(protConf001); 
        int [] numprotconf001 = countNumForwardAndDecoyHits(allpis);
        
        allpis.addAll(protConf01); 
        int [] numprotconf01 = countNumForwardAndDecoyHits(allpis);
        allpis.addAll(protConf02); 
        int [] numprotconf02 =  countNumForwardAndDecoyHits(allpis);
        allpis.addAll(protConf03); 
        int [] numprotconf03 =  countNumForwardAndDecoyHits(allpis);
        allpis.addAll(protConf05); 
        int [] numprotconf05 =  countNumForwardAndDecoyHits(allpis);
        allpis.addAll(protConf10); 
        int [] numprotconf10 =  countNumForwardAndDecoyHits(allpis);
  
/*
        int tspc0 = getTotalSpectrumCount(protConf0);
        int tspc001 = tspc0 + getTotalSpectrumCount(protConf001);
        int tspc01 = tspc001 + getTotalSpectrumCount(protConf01);
        int tspc02 = tspc01 + getTotalSpectrumCount(protConf02);
        int tspc03 = tspc02 + getTotalSpectrumCount(protConf03);
        int tspc05 = tspc03 + getTotalSpectrumCount(protConf05);
        int tspc10 = tspc05 + getTotalSpectrumCount(protConf10);
 
        int peptideseq0 = getNumPeptideSequences(protConf0);
        int peptideseq001 = peptideseq0 + getNumPeptideSequences(protConf001);
        int peptideseq01 = peptideseq001 + getNumPeptideSequences(protConf01);
        int peptideseq02 = peptideseq01 + getNumPeptideSequences(protConf02);
        int peptideseq03 = peptideseq02 + getNumPeptideSequences(protConf03);
        int peptideseq05 = peptideseq03 + getNumPeptideSequences(protConf05);
        int peptideseq10 = peptideseq05 + getNumPeptideSequences(protConf10);
*/
        // add more for peptide with no mods
        System.out.println("\n\n\nProtein Global FDR\t0\t0.001\t0.01\t0.02\t0.03\t0.05\t0.1");
        System.out.println("Number of Forward Peptide Hits\t" + numprotconf0[0] +  "\t" + numprotconf001[0] + "\t" + numprotconf01[0] + "\t" + numprotconf02[0] + "\t" + numprotconf03[0] + "\t" + numprotconf05[0] + "\t" + numprotconf10[0]);
        System.out.println("Number of Decoy Peptide Hits\t" + numprotconf0[1] +  "\t" + numprotconf001[1] + "\t" + numprotconf01[1] + "\t" + numprotconf02[1] + "\t" + numprotconf03[1] + "\t" + numprotconf05[1] + "\t" + numprotconf10[1]);
        //System.out.println("Number of Peptide Sequences\t" + peptideseq0 +  "\t" + peptideseq001 + "\t" + peptideseq01 + "\t" + peptideseq02 + "\t" + peptideseq03 + "\t" + peptideseq05 + "\t" + peptideseq10);
        System.out.println("Number of Forward Peptide Sequences\t" + numprotconf0[2] +  "\t" + numprotconf001[2] + "\t" + numprotconf01[2] + "\t" + numprotconf02[2] + "\t" + numprotconf03[2] + "\t" + numprotconf05[2] + "\t" + numprotconf10[2]);
        System.out.println("Number of Decoy Peptide Sequences\t" + numprotconf0[3] +  "\t" + numprotconf001[3] + "\t" + numprotconf01[3] + "\t" + numprotconf02[3] + "\t" + numprotconf03[3] + "\t" + numprotconf05[3] + "\t" + numprotconf10[3]);
        //System.out.println("Number of Spectrum Count\t" + tspc0 +  "\t" + tspc001 + "\t" + tspc01 + "\t" + tspc02 + "\t" + tspc03 + "\t" + tspc05 + "\t" + tspc10 + "\n\n\n");
        System.out.println("Number of Forward Hit Spectrum Count\t" + numprotconf0[4] +  "\t" + numprotconf001[4] + "\t" + numprotconf01[4] + "\t" + numprotconf02[4] + "\t" + numprotconf03[4] + "\t" + numprotconf05[4] + "\t" + numprotconf10[4]);
        System.out.println("Number of Decoy Hit Spectrum Count\t" + numprotconf0[5] +  "\t" + numprotconf001[5] + "\t" + numprotconf01[5] + "\t" + numprotconf02[5] + "\t" + numprotconf03[5] + "\t" + numprotconf05[5] + "\t" + numprotconf10[5]);


        System.out.println("\n\n\nCummulative Number of Protein Groups Identified");
        System.out.println(getCummulativeProteinCounts(allpgs));
        
        pepFile.close();
    // the allpgs should be sorted when this function is called, scoretype 2 for sumZScore, 5 for averageZScore, 6 for confidenceSum 
    }
    private static void calcIsoformScore(ArrayList<ProteinGroup> allpgs, int scoretype) {
       
        double [] localfdr = new double[allpgs.size()];
        int [] globlefdr = new int[allpgs.size()];
        boolean [] reversehits = new boolean[allpgs.size()];
        int numreversepgs = 0;
        int numpgs = 0;
//System.out.println("Number of Protein Groups Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tIsReverseHit");
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            numpgs++;
            if(pg.isReverseHit()) {
                reversehits[i] = true;
                globlefdr[i] = ++numreversepgs;
                
            }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
        }
        int numcloseby = 100; 
        for(int i = 0; i < reversehits.length; i++) {
            int temprev = 0;
            for(int j = i - numcloseby; j < i+numcloseby; j++) {
                if(j > -1 && j < reversehits.length && reversehits[j]) {
                    temprev++;
                }
            }
            localfdr[i] = temprev;
        }

       for(int i = 0; i < localfdr.length; i++) {
           int before = i >= numcloseby? numcloseby : i;
           int after = localfdr.length - i >= numcloseby? numcloseby : (localfdr.length - i);
//System.out.print(i + "\t" + localfdr[i] + "\t");
           localfdr[i] = (localfdr[i]+0.0)/(1 + before + after);
           if(i > 0) {
//System.out.println(localfdr[i] + "\t" + localfdr[i-1]); 
               if(localfdr[i] < localfdr[i-1]) {       
                   localfdr[i] = localfdr[i-1];
               }
           }
       }
      
//System.out.println("Protein Group ID\tNumber of Members\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tProtein Length\tIsReverseHit\tRepresentative Protein accession\tDefline\tAll accessions\tAll Deflines");

        int numforward = 0;
        int numreverse = 0; 
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            if(pg.isReverseHit()) {
                numreverse++; 
            } else {

                numforward++;
            }
    // the allpgs should be sorted when this function is called, scoretype 2 for sumZScore, 5 for averageZScore, 6 for confidenceSum 
    // 7 for confidenceProduct 
            switch(scoretype) {
                case 5 : pg.setAverageZScoreConfidence(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 2 : pg.setSumZScoreConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
                default : pg.setSumZScoreConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
            }


            String repacc = pg.getRepresentativeAcc();
            String repdefline = pg.getRepresentativeDefline();
            String accs = pg.getAccs();
            String deflines = pg.getDeflines();
             
//System.out.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + numforward + "\t" + numreverse + "\t" + (numreverse/(0.0 + numforward)) + "\t" + (1-localfdr[i]) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.getRepresentativeLength() + "\t" + pg.isReverseHit() + "\t" + repacc + "\t" + repdefline + "\t" + accs + "\t" + deflines);
        }
    }

    private static void calcProteinConfidenceByIdentifiedTrypticPeptides(ArrayList<ProteinGroup> allpgs) { 
        ArrayList<ProteinGroup> pgwith0trypticpeptides = new ArrayList();
        ArrayList<ProteinGroup> pgwith1trypticpeptides = new ArrayList();
        ArrayList<ProteinGroup> pgwithmoretrypticpeptides = new ArrayList();
        for(Iterator<ProteinGroup> it = allpgs.iterator(); it.hasNext();) {
            ProteinGroup pg = it.next();
            int numidedtyppeps = pg.getRepresentative().getNumIdentifiedTrypticPeptides(); 
            //int numyppeps = pg.getRepresentative().getNumTrypticPeptides(); 
            if(numidedtyppeps < 1) {
                pgwith0trypticpeptides.add(pg);
            } else if(numidedtyppeps == 1) {
                pgwith1trypticpeptides.add(pg);
            } else {
                pgwithmoretrypticpeptides.add(pg);
            }
        }
        Collections.sort(pgwith0trypticpeptides, new ProteinGroupComparator(6));
System.out.println("\n\n\nPeptide with no tryptic peptides identified");
        calcProteinConfidence(pgwith0trypticpeptides, 8, 100); // assign     
        Collections.sort(pgwith1trypticpeptides, new ProteinGroupComparator(8));
System.out.println("\n\n\nPeptide with one tryptic peptides identified");
        calcProteinConfidence(pgwith1trypticpeptides, 8, 100); // assign     
        Collections.sort(pgwithmoretrypticpeptides, new ProteinGroupComparator(8));
System.out.println("\n\n\nPeptide with more tryptic peptides identified");
        calcProteinConfidence(pgwithmoretrypticpeptides, 8, 100); // assign     
        System.out.println("Number of protein groups with no tryptic peptides:\t" + pgwith0trypticpeptides.size());
        System.out.println("Number of protein groups with 1 tryptic peptides:\t" + pgwith1trypticpeptides.size());
        System.out.println("Number of protein groups with more than 1 tryptic peptides:\t" + pgwithmoretrypticpeptides.size());
    }


    // the allpgs should be sorted when this function is called, scoretype 2 for sumZScore, 5 for averageZScore, 6 for confidenceSum 
    // 7 for confidenceProduct 
    private static void calcProteinConfidence(ArrayList<ProteinGroup> allpgs, int scoretype, int numcloseby) {
       
        double [] localfdr = new double[allpgs.size()];
        int [] globlefdr = new int[allpgs.size()];
        boolean [] reversehits = new boolean[allpgs.size()];
        int numreversepgs = 0;
        int numpgs = 0;
//System.out.println("Number of Protein Groups Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tIsReverseHit");
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            numpgs++;
            if(pg.isReverseHit()) {
                reversehits[i] = true;
                globlefdr[i] = ++numreversepgs;
            }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
        }
        //int numcloseby = 100; 
        for(int i = 0; i < reversehits.length; i++) {
            int temprev = 0;
            for(int j = i - numcloseby; j < i+numcloseby; j++) {
                if(j > -1 && j < reversehits.length && reversehits[j]) {
                    temprev++;
                }
            }
            localfdr[i] = temprev;
        }

       for(int i = 0; i < localfdr.length; i++) {
           int before = i >= numcloseby? numcloseby : i;
           int after = localfdr.length - i >= numcloseby? numcloseby : (localfdr.length - i);
//System.out.print(i + "\t" + localfdr[i] + "\t");
           localfdr[i] = (localfdr[i]+0.0)/(1 + before + after);
           if(i > 0) {
//System.out.println(localfdr[i] + "\t" + localfdr[i-1]); 
               if(localfdr[i] < localfdr[i-1]) {       
                   localfdr[i] = localfdr[i-1];
               }
           }
       }
      
//System.out.println("Protein Group ID\tNumber of Members\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tProtein Length\tIsReverseHit\tRepresentative Protein accession\tDefline\tAll accessions\tAll Deflines");

        int numforward = 0;
        int numreverse = 0; 
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            if(pg.isReverseHit()) {
                numreverse++; 
            } else {

                numforward++;
            }
    // the allpgs should be sorted when this function is called, scoretype 2 for sumZScore, 5 for averageZScore, 6 for confidenceSum 
    // 7 for confidenceProduct 
            switch(scoretype) {
                case 5 : pg.setAverageZScoreConfidence(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 2 : pg.setSumZScoreConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
                case 8 : pg.setTrypticConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
                default : pg.setSumZScoreConfidence(1-localfdr[i]); break; // set sumZScoreConfidence 
            }


            String repacc = pg.getRepresentativeAcc();
            String repdefline = pg.getRepresentativeDefline();
            String accs = pg.getAccs();
            String deflines = pg.getDeflines();
             
//System.out.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + numforward + "\t" + numreverse + "\t" + (numreverse/(0.0 + numforward)) + "\t" + (1-localfdr[i]) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.getRepresentativeLength() + "\t" + pg.isReverseHit() + "\t" + repacc + "\t" + repdefline + "\t" + accs + "\t" + deflines);
        }
    }
   
    // the allpgs should be sorted when this function is called 
    private static void calcProteinConfidence(ArrayList<ProteinGroup> allpgs) {
       
        double [] localfdr = new double[allpgs.size()];
        int [] globlefdr = new int[allpgs.size()];
        boolean [] reversehits = new boolean[allpgs.size()];
        int numreversepgs = 0;
        int numpgs = 0;
//System.out.println("Number of Protein Groups Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tIsReverseHit");
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            numpgs++;
            if(pg.isReverseHit()) {
                reversehits[i] = true;
                globlefdr[i] = ++numreversepgs;
                
            }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
        }
        int numcloseby = 100; 
        for(int i = 0; i < reversehits.length; i++) {
            int temprev = 0;
            for(int j = i - numcloseby; j < i+numcloseby; j++) {
                if(j > -1 && j < reversehits.length && reversehits[j]) {
                    temprev++;
                }
            }
            localfdr[i] = temprev;
        }

       for(int i = 0; i < localfdr.length; i++) {
           int before = i >= numcloseby? numcloseby : i;
           int after = localfdr.length - i >= numcloseby? numcloseby : (localfdr.length - i);
//System.out.print(i + "\t" + localfdr[i] + "\t");
           localfdr[i] = (localfdr[i]+0.0)/(1 + before + after);
           if(i > 0) {
//System.out.println(localfdr[i] + "\t" + localfdr[i-1]); 
               if(localfdr[i] < localfdr[i-1]) {       
                   localfdr[i] = localfdr[i-1];
               }
           }
       }
      
//System.out.println("Protein Group ID\tNumber of Members\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tProtein Length\tIsReverseHit\tRepresentative Protein accession\tDefline\tAll accessions\tAll Deflines");

        int numforward = 0;
        int numreverse = 0; 
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup pg = allpgs.get(i);
            if(pg.isReverseHit()) {
                numreverse++; 
            } else {

                numforward++;
            }
            String repacc = pg.getRepresentativeAcc();
            String repdefline = pg.getRepresentativeDefline();
            String accs = pg.getAccs();
            String deflines = pg.getDeflines();
//System.out.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + numforward + "\t" + numreverse + "\t" + (numreverse/(0.0 + numforward)) + "\t" + (1-localfdr[i]) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.getRepresentativeLength() + "\t" + pg.isReverseHit() + "\t" + repacc + "\t" + repdefline + "\t" + accs + "\t" + deflines);
        }
    }
    // put subset protein group in the super group, pgs are sorted by number of peptides
    private static ArrayList<ProteinGroup> subsetProteins(ArrayList<ProteinGroup> allpgs) {
System.out.println("Number of ProteinGroups before grouping subsets " + allpgs.size());
        ArrayList<ProteinGroup> pgs = new ArrayList(20000);
        int numpgs = 0;
        for(int i = 0; i < allpgs.size(); i++) {
            ProteinGroup current = allpgs.get(i);
            if(current != null) {
                pgs.add(current);
                numpgs++;
                allpgs.set(i, null);
                HashSet<PeptideItem> currentset = current.getPeptideItems();
                for(int j = i+1; j < allpgs.size(); j++) {
                
                    ProteinGroup next = allpgs.get(j);
                    if(next != null) {

                        HashSet<PeptideItem> nextset = next.getPeptideItems();
                        //if(currentset.size() == nextset.size() && currentset.containsAll(nextset)) {
                        if(currentset.containsAll(nextset)) {
                            current.addSubset(next);
                            allpgs.set(j, null);
                        } 

                    }
                }
            }
        }        

System.out.println("Final number of protein groups after grouping subsets " + numpgs);
        return pgs;
    }


    // remove pg from pgs if pg does not contain any unique peptide sequence
    private static ArrayList<ProteinGroup> removeNoUniquePepitdeProteins(ArrayList<ProteinGroup> pgs, ProteinGroup pg) {
        pgs.remove(pg);
        HashSet<PeptideItem> settobechecked = pg.getPeptideItems();
        HashSet<PeptideItem> allset = new HashSet(1000000);
        for(Iterator<ProteinGroup> it = pgs.iterator(); it.hasNext();) {
            allset.addAll(it.next().getPeptideItems());
        }
        if(!allset.containsAll(settobechecked)) {
            pgs.add(pg);
        } else {

//System.out.println("Protein with no unique peptide removed: " + pg.getAccessions() + "\tNumber of peptide: " + pg.getNumPeptideItems() + "\tNumber of spectral count: " + pg.getTotalSpectrumCount());
        }

        return pgs;
    }
    private static ArrayList<ProteinGroup> removeNoUniquePepitdeProteins(ArrayList<ProteinGroup> allpgs) {

System.out.println("Number of protein groups before removing no unique peptide proteins " + allpgs.size());
        ArrayList<ProteinGroup> temp = new ArrayList(allpgs.size());
        for(int i = allpgs.size() - 1; i > -1; i--) {
            ProteinGroup pg = allpgs.get(i);
            removeNoUniquePepitdeProteins(temp, pg);
        }
System.out.println("Number of protein groups after removing no unique peptide proteins " + temp.size());
        return temp; 
    }

    private static void assignPeptideConfidence(ArrayList<PeptideItem> pis, int scoretype) {
        
        Collections.sort(pis, new PeptideItemComparator(scoretype)); //sort by  
        double [] localfdr = new double[pis.size()];
        int [] globlefdr = new int[pis.size()];
        boolean [] reversehits = new boolean[pis.size()];
        int numreversepgs = 0;
        int numpgs = 0;
//System.out.println("Number of Protein Groups Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tIsReverseHit");
        for(int i = 0; i < pis.size(); i++) {
            PeptideItem pg = pis.get(i);
            numpgs++;
            if(pg.isReverseHit()) {
                reversehits[i] = true;
                globlefdr[i] = ++numreversepgs;
                
            }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
        }
        int numcloseby = 100; 
        for(int i = 0; i < reversehits.length; i++) {
            int temprev = 0;
            for(int j = i - numcloseby; j < i+numcloseby; j++) {
                if(j > -1 && j < reversehits.length && reversehits[j]) {
                    temprev++;
                }
            }
            localfdr[i] = temprev;
        }

       for(int i = 0; i < localfdr.length; i++) {
           int before = i >= numcloseby? numcloseby : i;
           int after = localfdr.length - i >= numcloseby? numcloseby : (localfdr.length - i);
//System.out.print(i + "\t" + localfdr[i] + "\t");
           localfdr[i] = (localfdr[i]+0.0)/(1 + before + after);
           if(i > 0) {
//System.out.println(localfdr[i] + "\t" + localfdr[i-1]); 
               if(localfdr[i] < localfdr[i-1]) {       
                   localfdr[i] = localfdr[i-1];
               }
           }
       }
      
System.out.println("Protein Group ID\tNumber of Members\tNumber of Forward Protein Hit\tNumber of Reverse Protein Hit\tGloble False Positive Rate\tConfidence(1 - Local False Positive Rate)\tNumber of PeptideItem\tTotal Spectrum Count\tsumZScore\tProtein Length\tIsReverseHit\tRepresentative Protein accession\tDefline\tAll accessions\tAll Deflines");

        int numforward = 0;
        int numreverse = 0; 
        for(int i = 0; i < pis.size(); i++) {
            PeptideItem pi = pis.get(i);
            if(pi.isReverseHit()) {
                numreverse++; 
            } else {

                numforward++;
            }
            switch(scoretype) {
                case 1 : pi.addConfidenceScore(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 2 : pi.addConfidenceScore(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 3 : pi.addConfidenceScore(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 4 : pi.addConfidenceScore(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 5 : pi.addConfidenceScore(1-localfdr[i]); break; // set averageZScoreConfidence 
                case 8 : pi.addConfidenceScore(1-localfdr[i]); break; // set averageZScoreConfidence 
                default : pi.setConfidence(1-localfdr[i]);

//System.out.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + numforward + "\t" + numreverse + "\t" + (numreverse/(0.0 + numforward)) + "\t" + (1-localfdr[i]) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.getRepresentativeLength() + "\t" + pg.isReverseHit() + "\t" + repacc + "\t" + repdefline + "\t" + accs + "\t" + deflines);
                          break; // set sumZScoreConfidence 
 
            }


             
//System.out.println((i+1) + "\t" + pg.getNumProteinItems() + "\t" + numforward + "\t" + numreverse + "\t" + (numreverse/(0.0 + numforward)) + "\t" + (1-localfdr[i]) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getTotalSpectrumCount() + "\t" + pg.getSumZScore() + "\t" + pg.getRepresentativeLength() + "\t" + pg.isReverseHit() + "\t" + repacc + "\t" + repdefline + "\t" + accs + "\t" + deflines);
        }


    }

    private static void calcPeptideConfidenceScore(ArrayList<PeptideItem> pis) {
        System.out.println("Sorting peptides by best ZScore");
        //assignPeptideConfidence(pis, 1); // with best ZScore 
        assignPeptideConfidence(pis, 2); // with sumZSorr 
        //System.out.println("Sorting peptides by sumXCorr");
        //assignPeptideConfidence(pis, 3); // with sumXCorr 
        //assignPeptideConfidence(pis, 4); // with occurence 
        //assignPeptideConfidence(pis, 5); // with bestXCorr // does not help much 
//System.out.println("now by avgZScore");
        assignPeptideConfidence(pis, 8); // with average ZScore

        System.out.println("Sorting peptides by sumConfidence");
        assignPeptideConfidence(pis, 6); // with sumConfidence 

    }

     
    private static ArrayList<PeptideItem> filterPeptideItems(HashMap<String, PeptideItem> seqcharge2peptideitem, double peptidefdr) {

        HashSet<String> uniquepeptideseqs = new HashSet(1000000);
        HashSet<String> uniquereversedpeptideseqs = new HashSet(1000000);
        HashSet<String> uniqueforwardpeptideseqs = new HashSet(1000000);

        HashSet<String> filtereduniquereversedpeptideseqs = new HashSet(1000000);
        HashSet<String> filtereduniqueforwardpeptideseqs = new HashSet(1000000);

        ArrayList<PeptideItem> filteredPeptideItems = new ArrayList();
        ArrayList<PeptideItem> pis = new ArrayList();
 
        pis.addAll(seqcharge2peptideitem.values());

        // now using confidence for peptides
        calcPeptideConfidenceScore(pis); 
        Collections.sort(pis, new PeptideItemComparator(6));
        //Collections.sort(pis);

        int numforward = 0;
        int numreverse = 0;
        System.out.println("Number of Forward Hits\tNumber of Reverse Hit\tFalse Positive Rate\tSumZScore\tOccurrence\tTotal Spectrum Count\tNumber of Unique Peptide Sequence added\tSequence\tCharge\tIsReverseHit");
        for(Iterator<PeptideItem> it = pis.iterator(); it.hasNext();) {
            PeptideItem pi = it.next();
            
            uniquepeptideseqs.add(pi.getSequence());
            if(pi.isReverseHit()) {
                numreverse++;
                uniquereversedpeptideseqs.add(pi.getSequence());
            } else {
                numforward++;
                uniqueforwardpeptideseqs.add(pi.getSequence()); 
            }

            double fdr = uniquereversedpeptideseqs.size()/(uniqueforwardpeptideseqs.size()+0.0);

            if(fdr <= peptidefdr) {
                filteredPeptideItems.add(pi);
                if(pi.isReverseHit()) {
                    filtereduniquereversedpeptideseqs.add(pi.getSequence());
                } else {
                    filtereduniqueforwardpeptideseqs.add(pi.getSequence()); 
                }
                
            }
            //System.out.println(numforward + "\t" + numreverse + "\t" + (numreverse+0.0)/(numforward) + "\t" + pi.getSumZScore() + "\t" + pi.getOccurrence() + "\t" + pi.getTotalSpectrumCount() + "\t" + uniquepeptideseqs.size() + "\t" + pi.getSequence() + "\t" + pi.getCharge() + "\t" + pi.isReverseHit());
        }
       
        System.out.println("\n\nTotal Number of PeptideItems Accepted: " + filteredPeptideItems.size() + "\tout of  Total Number of PeptideItems of " + (numreverse + numforward)+"\tacceptedUniqueReversedSeq: " + filtereduniquereversedpeptideseqs.size() + "\tacceptedUniqueForwardSeq: " + filtereduniqueforwardpeptideseqs.size()); 
        return filteredPeptideItems; 
    }

    // adjust protein score based on the average score of reversed ids with similar protein length
    public static void assignProteinScore(HashSet<ProteinItem> proteins) {

        HashMap<ProteinItem, Double> proteinitem2avgreversescore = new HashMap(1000000); // average of reversed protein scores
        ArrayList<ProteinItem> prots = new ArrayList(proteins);
        Collections.sort(prots, new ProteinItemComparator(1)); //sort by protein length
        //Collections.reverse(prots); // short to length
        int windowsize = 200;
        ArrayList<ProteinItem> templist = new ArrayList(windowsize);
        int numprocessed = 0;
        double sumreversedscore = 0; 
        int numreversed = 0;
        int numlengthbin = 100;
        int [] lengthfreq = new int[numlengthbin];
        int [] reversedlengthfreq = new int[numlengthbin];
        int totalforwardprot = 0;       
        int totalreversedprot = 0;       
 
        for(Iterator<ProteinItem> it = prots.iterator(); it.hasNext();) {
            ProteinItem pi = it.next();
            int protlength = pi.getFasta().getLength();
            if(pi.isReverseHit()) { 
                totalreversedprot++;
                if(protlength >= 10000) { 
                    reversedlengthfreq[numlengthbin-1]++;
                } else {
                    reversedlengthfreq[protlength/100]++;
                }
            } else {
                totalforwardprot++;
                if(protlength >= 10000) { 
                    lengthfreq[numlengthbin-1]++;
                } else {
                    lengthfreq[protlength/100]++;
                }

            }
            
            // correct by the protein score by protein length
System.out.print("protein length: " + pi.getFasta().getLength());
            String ac = pi.getFasta().getLongAccession();
            double score = pi.getSumZScore();
            numprocessed++;
            if(ac.startsWith("Revers")) {
                sumreversedscore +=  pi.getSumZScore();
                numreversed++;
            } 
            templist.add(pi);
            if(templist.size() == windowsize || numprocessed == prots.size() ) {
                double avgreverscore = sumreversedscore/numreversed;
System.out.println("avgreverscore: " + avgreverscore + "\tnumprocessed: " + numprocessed);
                Iterator<ProteinItem> it1 = templist.iterator();
                while(it1.hasNext()) {
                    ProteinItem temppi = it1.next();
                    //temppi.setProteinScore(temppi.getSumZScore() - avgreverscore);
                    temppi.setProteinScore(temppi.getBestPeptideScore());
  
                }

                templist = new ArrayList(windowsize);
                sumreversedscore = 0;
                numreversed = 0;

            }
          
        }
        StringBuffer lensb = new StringBuffer();
        StringBuffer freqsb = new StringBuffer();
        StringBuffer reversedfreqsb = new StringBuffer();
        StringBuffer relativefreqsb = new StringBuffer();
        StringBuffer relativereversedfreqsb = new StringBuffer();
        lensb.append("Protein Length");
        freqsb.append("Forward Protein Frequency");
        reversedfreqsb.append("Reversed Protein Frequency");
        relativefreqsb.append("Forward Protein Relative Frequency");
        relativereversedfreqsb.append("Reversed Protein Relative Frequency");
        
        for(int i = 0; i < numlengthbin; i++) {
            lensb.append("\t" + i*100 + "-" + (i+1)*100);
            freqsb.append("\t" + lengthfreq[i]);
            relativefreqsb.append("\t" + lengthfreq[i]/(0.0+totalforwardprot));
            reversedfreqsb.append("\t" + reversedlengthfreq[i]);
            relativereversedfreqsb.append("\t" + reversedlengthfreq[i]/(0.0+totalreversedprot));
            
        } 
        System.out.println(lensb);
        System.out.println(freqsb);
        System.out.println(reversedfreqsb);
        System.out.println(relativefreqsb);
        System.out.println(relativereversedfreqsb);
    }
/*
    // get the conterpart protein - ie the same protein, reversed and forward copy
    public static void assignProteinScore(HashSet<ProteinItem> proteins) {
        HashMap<String, ProteinItem> ac2protein = new HashMap(1000000);
        HashMap<String, ProteinItem> reverseac2protein = new HashMap(1000000);
       
        for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
            ProteinItem pi = it.next();
            String ac = pi.getFasta().getLongAccession();
            
            ac2protein.put(ac, pi);
          
        }

        for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
            ProteinItem pi = it.next();
            String ac = pi.getFasta().getLongAccession();
            
            ProteinItem counterpart = null;
            String counterpartac = null;
            if(ac.startsWith("Revers")) { 
                counterpartac = ac.split("Reverse_")[1];
            } else {
                
                counterpartac = ("Reverse_") + ac;
            }

            counterpart = ac2protein.get(counterpartac);
            double score = pi.getSumZScore();
            if(counterpart != null) {
                score -= counterpart.getSumZScore();
System.out.println("Counterpart found for " + ac + "\t counterpartac: " + counterpartac );
            } else {
System.out.println("Counterpart not found for " + ac + "\t counterpartac: " + counterpartac );
            }
            pi.setProteinScore(score);
        }
        
    }
*/

    // in case the pepseq is found in both forward and reverse protein, return only the forward protein
    public static ArrayList<Fasta> peptideseq2Fastas(String pepseq, PrefixDb pdb) {
        ArrayList<Fasta> ffastas = new ArrayList(); // for forward hits
        ArrayList<Fasta> rfastas = new ArrayList(); // for reversed hits
         
        Iterator<Fasta> prots = pdb.peptideseq2Fastas(pepseq).iterator();
        while(prots.hasNext()) {
            Fasta f = prots.next();
            String  protseq = f.getSequence();
            //String myac = f.getAccession();
            
            if(protseq.indexOf(pepseq) != -1) {
                if(f.isReversed()) {
                    rfastas.add(f);    
                } else {
                    ffastas.add(f);
                }
            }
        }
        if(ffastas.size() == 0) {
            return rfastas;
        } else {
            return ffastas;
        }
    }
    public static HashSet<ProteinItem> mapPeptide2Proteins(String database,  ArrayList<PeptideItem> acceptedpeptides) throws IOException {
        HashSet<ProteinItem> proteins = new HashSet(1000000);
        HashMap<String, ProteinItem> ac2protitem = new HashMap(1000000);

        PrefixDb pdb = new PrefixDb(database);

        int numpeptideprocessed = 0;
        for(Iterator<PeptideItem> it = acceptedpeptides.iterator(); it.hasNext();) {
            PeptideItem pi = it.next();
            //String pepseq = pi.getSequence();
            String pepseq = pi.getSeqWithoutMod();
            //Iterator<Fasta> prots = pdb.peptideseq2Fastas(pepseq).iterator();
            Iterator<Fasta> prots = peptideseq2Fastas(pepseq, pdb).iterator();
            while(prots.hasNext()) {
                Fasta f = prots.next();
                String  protseq = f.getSequence();
                String myac = f.getAccession();
                ProteinItem proti = ac2protitem.get(myac);
                if(protseq.indexOf(pepseq) != -1) {
                    if(proti == null) {
                        proti = new ProteinItem(f);
                        proteins.add(proti);
                        ac2protitem.put(myac, proti);

                    }

                    pi.addProteinItem(proti);
                    proti.addPeptideItem(pi);
                    
                }
                
            }

//System.out.print("Number of Peptides Processed: " + ++numpeptideprocessed + "\r");
        }
        
        int totalnumpeptides = 0;
        for(Iterator<ProteinItem> it = proteins.iterator(); it.hasNext();) {
            ProteinItem pi = it.next();
            totalnumpeptides += pi.getNumPeptides();

        }
        int numproteinadded = proteins.size();
System.out.println("\nNumber of protein items added " + proteins.size() + "\tAverage number of peptide per protein: " + totalnumpeptides/(0.0 + numproteinadded) );
for(Iterator<PeptideItem> it = acceptedpeptides.iterator(); it.hasNext();) {
    PeptideItem pi = it.next();
    if(pi.isReverseHit() && pi.containsBothForwardAndReversedHits()) {
System.out.println("Peptide with found in both forward and reversed: " + pi.getSequence() + "\tsumXCorr: " + pi.getSumXCorr() + "\tbest ZScore: " + pi.getBestZScore() + "\tOccurrence: " + pi.getOccurrence() + "\tsumZScore: " + pi.getSumZScore());
    }
}
/*


        FileInputStream fis = new FileInputStream(new File(database));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        int count = 0;
        int numproteinadded = 0;
        int totalnumpeptides = 0;
        while(fastas.hasNext()) {
           
            Fasta f = fastas.next();
            count++;
            String protseq = f.getSequence();
            ProteinItem proti = null;
            String myac = f.getAccession();
            Iterator<PeptideItem> it = acceptedpeptides.iterator();
System.out.print("processing " + myac + " " + count + "\r");
            int numpeptides = 0;
            while(it.hasNext()) {
                PeptideItem pi = it.next();
                //String pepseq = pi.getSequence();
                String pepseq = pi.getSeqWithoutMod();
//System.out.println("processing peptide sequence " + pepseq);
                if(protseq.indexOf(pepseq) != -1) {
                    numpeptides++;
                    if(proti == null) {
                        proti = new ProteinItem(f);
                        proteins.add(proti);
                        numproteinadded++;

                    }

                    pi.addProteinItem(proti);
                    proti.addPeptideItem(pi);
                }

            }
           
            totalnumpeptides += numpeptides; 
        }
*/
        return proteins;
    }
    
    private static String getDatabase()
    {
        String databaseName = "";
        try {
            BufferedReader br = new BufferedReader(new FileReader(folders.get(0)+File.separator+"sequest.params"));
            while(true)
            {
                String currentLine = br.readLine();
                if(currentLine.startsWith("database_name"))
                {
                    databaseName = currentLine.split("=")[1];
                    break;
                }
                
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ProteinInference.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ProteinInference.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        return databaseName;
    }
    public static void main(String[] args) throws Exception{

       
        double peptidefdr = 0.05; 
//       System.out.println(args[0]);

        peptideScoreType = 1;
	BufferedReader br = null;
        
        //if(args.length < 2) {
        //    help();
        //}
        //ac2Fasta =  getProteinMap(database);
        int copacount = 0;
	try {
            if(args.length < 2) {
                System.out.println(HELP);
                System.exit(0);

            }
            if(args[0].contains(","))
            {
                if(args.length > 1) {
                    database = args[1];
                }
                if(args.length > 2) {
                    peptidefdr = Double.parseDouble(args[2]);
                }
                if(args.length > 0) {
                    for(String fileName :args[0].split(",") )
                    {
                        folders.add(fileName);
  //                      System.out.println(fileName);
                    }
                }
                
            
            
            }
            else
            {

                if(args.length > 0) {
                    dtaselectfiles = args[0];
                }
                if(args.length > 1) {
                    database = args[1];
                }
                if(args.length > 2) {
                    peptidefdr = Double.parseDouble(args[2]);
                }

                br = new BufferedReader(new FileReader(dtaselectfiles));
                String line = null;
                line = br.readLine();
                while(line != null) {
                    line = line.trim();
                    if(!"".equals(line)) {
                        File dtafile = new File(line);
                        if(dtafile.isFile() || dtafile.isDirectory()) {
                            folders.add(line); 
                        } else {

                            System.out.println(dtafile + " does not exist");
                        }

                    } 
                    line = br.readLine();
                } 
            br.close();
            }
            
            if(args.length == 4)
                outputPath = args[3];
            if(database.equals("."))
            {
                database = getDatabase();
            }
//            System.out.println("---------------------------"+outputPath);

            int numpeptide = 0;
            HashSet<String> uniquepepitdes = new HashSet(10000000);
            HashSet<String> uniquepepitdesignorecharge = new HashSet(10000000);
            HashMap<String, PeptideItem> seqcharge2peptideitem = new HashMap(1000000);
//ps.println("Number of experiment: " + folders.size());
            for(int i = 0; i < folders.size(); i++) {
                String folder = folders.get(i);
                if(!"".equals(folder)) {
                    System.out.println("Now processing " + folder);
                    String dtafilter = folder;
                    if(!dtafilter.endsWith("DTASelect-filter.txt")) {
                        dtafilter = folder + "/DTASelect-filter.txt";
                    }


                    DTASelectFilterReader reader = new DTASelectFilterReader(dtafilter);
                    HashSet<String> uniquepeptides = new HashSet(1000000);
                    for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
                        Protein p = itr.next();
                        String proteinid = p.getAccession();
                        boolean isreversehit = proteinid.startsWith("Reverse");
//System.out.println("protien id: " + proteinid);
                        
                        for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                            numpeptide++;
                            Peptide peptide = pepItr.next();
                            //String pepseq = peptide.getSeqWithNoModification();
                            //String seq = peptide.getSequence();
                            String seq = peptide.getMidSeq();
                            String filename = peptide.getFileName();
                            String scan = peptide.getLoScan();
                            String charge = peptide.getChargeState();
                            int chargestate = Integer.parseInt(charge);
                            int scanInt = Integer.parseInt(scan);
                            
                            double spvalue = peptide.getSpScoreValue();
                            String seqcharge = seq + charge; 
                            uniquepepitdes.add(seqcharge);
                            uniquepepitdesignorecharge.add(seq);
                            PeptideItem pi = seqcharge2peptideitem.get(seqcharge);
                            if(pi == null) {
                                pi = new PeptideItem(seq, chargestate, isreversehit, peptideScoreType, folders.size());
                                seqcharge2peptideitem.put(seqcharge, pi);
                            }
//                            else
//                            {
//                                PeptideItem newPI =  new PeptideItem(seq, chargestate, isreversehit, peptideScoreType, folders.size());
//                                int newSPC = newPI.get + pi.getTotalSpectrumCount();
//                                pi.setTotalSpectrumCount(newSPC);
//                                seqcharge2peptideitem.put(seqcharge, pi);
//                            }
                            
                            if(!uniquepeptides.contains(seqcharge)) {
                                pi.addPeptide(peptide, i); // to ensure the seqcharge only be added only for each DTASelect-filter.txt file
                                uniquepeptides.add(seqcharge);
                            }
                            
                            //System.out.println("File name: " + filename + "\t" + scan);

                        }
                    }
                    System.out.println("finished reading DTASelect-filter.txt");
                    reader.close();
                } 
            } 

            System.out.println("Number of Peptide Identified: " + numpeptide +  "\tNumber of unique peptides: " + uniquepepitdes.size() + "\tNumber of unique peptides ignore charge: " + uniquepepitdesignorecharge.size());
            System.out.println("Starting to map peptides to proteins");

            long starttime = System.nanoTime(); 
            ArrayList<PeptideItem> acceptedPeptides = filterPeptideItems(seqcharge2peptideitem, peptidefdr);
            long endtime = System.nanoTime();
            System.out.println("Time used to filter peptide items: " + (endtime - starttime)/1000000);
            starttime = endtime; 

            //HashSet<ProteinItem> protitems = mapPeptide2Proteins(database, seqcharge2peptideitem);
            HashSet<ProteinItem> protitems = mapPeptide2Proteins(database, acceptedPeptides);
            endtime = System.nanoTime();
            System.out.println("Time used to map peptide to proteins: " + (endtime - starttime)/1000000);
            starttime = endtime; 

            // porteinScores will not be used in the final confidence, so this statement might be deleted
            //assignProteinScore(protitems);
             
            System.out.println("Finished mapping peptides to proteins");
            System.out.println("Start to group proteins");

            ArrayList<ProteinGroup> pgs = groupProteins(protitems);
            System.out.println("Finished grouping proteins");

            endtime = System.nanoTime();
            System.out.println("Time used to group proteins: " + (endtime - starttime)/1000000);
            starttime = endtime; 

            // to properly remove subset proteins, it is important for sort protein groups by pepitde number
            Collections.sort(pgs, new ProteinGroupComparator(4)); // sort by peptide number

            ArrayList<ProteinGroup> finalpgs = subsetProteins(pgs);

            endtime = System.nanoTime();
            System.out.println("Time used to subset proteins: " + (endtime - starttime)/1000000);
            starttime = endtime; 

            Collections.sort(finalpgs); 
            System.out.println("Before removing proteins with no unique pepetides: " + finalpgs.size());
            ArrayList<ProteinGroup> pgsafterremovingnouniquepepitde =  removeNoUniquePepitdeProteins(finalpgs);
            System.out.println("After removing proteins with no unique pepetides: " + pgsafterremovingnouniquepepitde.size());

            endtime = System.nanoTime();
            System.out.println("Time used to remove protein groups with no unique peptides: " + (endtime - starttime)/1000000);
            starttime = endtime; 
        
            assignPeptideItemDecoyStatus(pgsafterremovingnouniquepepitde);

            Collections.sort(pgsafterremovingnouniquepepitde);
            calcProteinConfidence(pgsafterremovingnouniquepepitde);


            System.out.println("\n\n\nSorting protein by sumZScore and calculating protein confidence based on sumZScore");
            Collections.sort(pgsafterremovingnouniquepepitde, new ProteinGroupComparator(2));
            calcProteinConfidence(pgsafterremovingnouniquepepitde, 2, 100);


            System.out.println("\n\n\nSorting protein by AverageZScore and calculating protein confidence basedon averageZScore");
            Collections.sort(pgsafterremovingnouniquepepitde, new ProteinGroupComparator(5));
            calcProteinConfidence(pgsafterremovingnouniquepepitde, 5, 100);
          
          
            //System.out.println("\n\n\nSorting protein by trypitc peptide fraction");
            //does not have positive impact, so disabled
            //calcProteinConfidenceByIdentifiedTrypticPeptides(pgsafterremovingnouniquepepitde);

            System.out.println("\n\n\nSorting protein by confidenceSum"); // this works the best
            Collections.sort(pgsafterremovingnouniquepepitde, new ProteinGroupComparator(6));
            //calcProteinConfidence(pgsafterremovingnouniquepepitde, 6);
            outputProteinConfidence(pgsafterremovingnouniquepepitde, 6);

            //System.out.println("Sorting protein by confidenceProduct");
            //Collections.sort(pgsafterremovingnouniquepepitde, new ProteinGroupComparator(7));
            //calcProteinConfidence(pgsafterremovingnouniquepepitde, 7);
          
            outputPeptideWithProteinConfidence(pgsafterremovingnouniquepepitde);
            endtime = System.nanoTime();
            System.out.println("Time used to calculate protein confidence: " + (endtime - starttime)/1000000);
            starttime = endtime; 

            //ps.close();
            //System.out.println("Number of Peptide Identified: " + numpeptide + "\tNumber of unique peptides: " + uniquepepitdes.size() + "\tNumber from Map: " + seqcharge2BestPeptide.size());
            

	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }

     
    private static String getUniProtAccession(String longac) {
        String [] arr = longac.split("\\|");
        String ac = arr[0];
        if(arr.length > 1) {
            if(longac.startsWith("Reverse")) {
                ac = arr[0] + "|" + arr[1];
            } else {
                ac = arr[1];
            }
        }
        return ac;
    }
}

