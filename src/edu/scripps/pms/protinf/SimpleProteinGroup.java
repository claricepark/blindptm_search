package edu.scripps.pms.protinf;

import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;

import edu.scripps.pms.util.seq.PrefixDb;


// 
public class SimpleProteinGroup implements Comparable {

    private ArrayList<SimpleProteinItem> proteins = new ArrayList();
    private ArrayList<SimpleProteinGroup> subset = new ArrayList();
    private SimpleProteinItem representative = null;
    

    //private ArrayList<Peptide> [] peptides = null;
    private boolean isAssignedToPeptideItem = false;
    private boolean isNoUniquePeptideGroup = false;
   

    public SimpleProteinItem getRepresentative() {
        return representative;
    } 

    public int compareTo(Object o) {
        //int num1 = peptides.size();
        //int num2 = ((SimpleProteinItem)o).getPeptideItems().size();
        //double num1 = getRepresentative().getSumZScore();
        //double num2 = ((SimpleProteinGroup)o).getRepresentative().getSumZScore();
        double num1 = getRepresentative().getNumPeptides();
        double num2 = ((SimpleProteinGroup)o).getRepresentative().getNumPeptides();
        if(num1 == num2) {
            return 0;
        } else if(num1 < num2) {
            return 1;
        } else {
            return -1;
        }

    }

    public void addProteinItem(SimpleProteinItem pit) {
        // check representative
        if(representative == null) {
            representative = pit;
        } else {
            if(pit.getLength() < representative.getLength()) {
                representative = pit;
            }
        }
        proteins.add(pit);
    }
    public static HashSet<SimpleProteinItem> mapPeptide2Proteins(String database,  HashSet<String> acceptedpeptides) throws IOException {
        HashSet<SimpleProteinItem> proteins = new HashSet(1000000);
        HashMap<String, SimpleProteinItem> ac2protitem = new HashMap(1000000);

        PrefixDb pdb = new PrefixDb(database);

        int numpeptideprocessed = 0;
        for(Iterator<String> it = acceptedpeptides.iterator(); it.hasNext();) {
            //String pepseq = pi.getSequence();
            String pepseq = it.next();
            //Iterator<Fasta> prots = pdb.peptideseq2Fastas(pepseq).iterator();
            Iterator<Fasta> prots = ProteinInference.peptideseq2Fastas(pepseq, pdb).iterator();
            while(prots.hasNext()) {
                Fasta f = prots.next();
                String  protseq = f.getSequence();
                String myac = f.getAccession();
                SimpleProteinItem proti = ac2protitem.get(myac);
                if(protseq.indexOf(pepseq) != -1) {
                    if(proti == null) {
                        proti = new SimpleProteinItem(f);
                        proti.addPeptide(pepseq);
                        proteins.add(proti);
                        ac2protitem.put(myac, proti);

                    }

                }
                
            }
           
//System.out.print("Number of Peptides Processed: " + ++numpeptideprocessed + "\r");
        }
        
        int totalnumpeptides = 0;
        for(Iterator<SimpleProteinItem> it = proteins.iterator(); it.hasNext();) {
            SimpleProteinItem pi = it.next();
            totalnumpeptides += pi.getNumPeptides();

        }
        int numproteinadded = proteins.size();
System.out.println("\nNumber of protein items added " + proteins.size() + "\tAverage number of peptide per protein: " + totalnumpeptides/(0.0 + numproteinadded) );
        return proteins;
    }

    public boolean isReverseHit() {
        for(Iterator<SimpleProteinItem> it = proteins.iterator(); it.hasNext();) {
            if(!it.next().isReverseHit()) return false;
        }
        return true;
    }
    
    public static ArrayList<SimpleProteinGroup> groupProteins(HashSet<SimpleProteinItem>  pis) {
//System.out.println("Number of ProteinItems before grouping " + pis.size());
        ArrayList<SimpleProteinGroup> pgs = new ArrayList(100000);
        ArrayList<SimpleProteinItem> temppis = new ArrayList();
        temppis.addAll(pis);
        Collections.sort(temppis);
        int numpgs = 0;
        int numreversepgs = 0;
//System.out.println("GroupProtein Added\tNumReverseHit\tFDR\tNumber of PeptideItem\tsumZScore\tIsReverseHit");
        for(int i = 0; i < temppis.size(); i++){
            SimpleProteinItem current = temppis.get(i);
            if(current != null) {
                SimpleProteinGroup pg = new SimpleProteinGroup();
                numpgs++;
                pg.addProteinItem(current);
                if(pg.isReverseHit()) {
                    numreversepgs++;
                }
//System.out.println(numpgs + "\t" + numreversepgs + "\t" + numreversepgs/(0.0 + numpgs) + "\t" + pg.getNumPeptideItems() + "\t" + pg.getSumZScore() + "\t" + pg.isReverseHit());
                pgs.add(pg);
                temppis.set(i, null);
                HashSet<String> currentset = current.getPeptideItems();
                for(int j = i+1; j < temppis.size(); j++) {
                
                    SimpleProteinItem next = temppis.get(j);
                    if(next != null) {

                        HashSet<String> nextset = next.getPeptideItems();
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
    public void addSubset(SimpleProteinGroup pg) {
        subset.add(pg);
    }
    public HashSet<String> getPeptideItems() {
        return getRepresentative().getPeptideItems(); 
    }

    // put subset protein group in the super group, pgs are sorted by number of peptides
    private static ArrayList<SimpleProteinGroup> subsetProteins(ArrayList<SimpleProteinGroup> allpgs) {
System.out.println("Number of SimpleProteinGroups before grouping subsets " + allpgs.size());
        ArrayList<SimpleProteinGroup> pgs = new ArrayList(20000);
        int numpgs = 0;
        for(int i = 0; i < allpgs.size(); i++) {
            SimpleProteinGroup current = allpgs.get(i);
            if(current != null) {
                pgs.add(current);
                numpgs++;
                allpgs.set(i, null);
                HashSet<String> currentset = current.getPeptideItems();
                for(int j = i+1; j < allpgs.size(); j++) {
                
                    SimpleProteinGroup next = allpgs.get(j);
                    if(next != null) {

                        HashSet<String> nextset = next.getPeptideItems();
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
    private static ArrayList<SimpleProteinGroup> labelNoUniquePepitdeProteins(ArrayList<SimpleProteinGroup> allpgs) {

System.out.println("Number of protein groups before removing no unique peptide proteins " + allpgs.size());
        //ArrayList<SimpleProteinGroup> temp = new ArrayList(allpgs.size());
        ArrayList<SimpleProteinGroup> temp = allpgs;
        for(int i = allpgs.size() - 1; i > -1; i--) {
            SimpleProteinGroup pg = allpgs.get(i);
            labelNoUniquePepitdeProteins(temp, pg);
        }
System.out.println("Number of protein groups after removing no unique peptide proteins " + temp.size());
        return temp; 
    }
    private static ArrayList<SimpleProteinGroup> labelNoUniquePepitdeProteins(ArrayList<SimpleProteinGroup> pgs, SimpleProteinGroup pg) {
        
        HashSet<String> settobechecked = pg.getPeptideItems();
        HashSet<String> allset = new HashSet(1000000);
        for(Iterator<SimpleProteinGroup> it = pgs.iterator(); it.hasNext();) {
            SimpleProteinGroup temppg = it.next();
            if(temppg != pg && !temppg.isNoUniquePeptideGroup()) { 
                allset.addAll(temppg.getPeptideItems());
            }
        }
        if(!allset.containsAll(settobechecked)) {
            
            pg.isNoUniquePeptideGroup(false);
        } else {
            pg.isNoUniquePeptideGroup(true);
//System.out.println("Protein with no unique peptide removed: " + pg.getAccessions() + "\tNumber of peptide: " + pg.getNumPeptideItems() + "\tNumber of spectral count: " + pg.getTotalSpectrumCount());
        }

        return pgs;
    }
    public boolean isNoUniquePeptideGroup() {
        return isNoUniquePeptideGroup;
    }

    public void isNoUniquePeptideGroup(boolean nouniquepepites) {
        isNoUniquePeptideGroup = nouniquepepites;
    }

    public static ArrayList<SimpleProteinGroup> pepseqs2SimpleProteinGroups(HashSet<String> peptides, String fastafile) throws IOException {
          
        HashSet<SimpleProteinItem> protitems = mapPeptide2Proteins(fastafile, peptides);
        ArrayList<SimpleProteinGroup> pgs = groupProteins(protitems);
System.out.println("Number of pgs before subset: " + pgs.size());
        // to properly remove subset proteins, it is important for sort protein groups by pepitde number
        Collections.sort(pgs); 
        ArrayList<SimpleProteinGroup> finalpgs = subsetProteins(pgs);
        Collections.sort(finalpgs); 
        System.out.println("Before removing proteins with no unique pepetides: " + finalpgs.size());
        ArrayList<SimpleProteinGroup> pgsafterremovingnouniquepepitde =  labelNoUniquePepitdeProteins(finalpgs);

        return pgsafterremovingnouniquepepitde;
    }
    public static void main(String [] args) throws IOException {
        HashSet<String> peptides = new HashSet();
        peptides.add("PDVLTTGGGN");
        peptides.add("SLTVGPRG");
        peptides.add("RDDIA");
        peptides.add("ETWGINHVNEDGTIEI");
        peptides.add("SVGEKDWQER");


   
        ArrayList<SimpleProteinGroup> spg = pepseqs2SimpleProteinGroups(peptides, "/lustre/people/applications/yates/dbase/17PM_pombe030305.fasta");
    }
}

