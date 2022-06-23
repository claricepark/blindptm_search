package edu.scripps.pms.protinf;

import java.util.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.util.dtaselect.Peptide;



// 
public class SimpleProteinItem implements Comparable {
    private HashSet<String> peptides = new HashSet();
    private Fasta fasta = null;
   
    public SimpleProteinItem(Fasta f) {
        fasta = f;
    }


    public void addPeptide(String p) {
        peptides.add(p);
    }
    public HashSet<String> getPeptideItems() {
        return peptides;
    }
    public int getLength() {
        return fasta.getLength();
    }


    public boolean isReverseHit() {
        return fasta.getAccession().startsWith("Revers");
    }
    public int getNumPeptides() {
        return peptides.size();
    }
    public int compareTo(Object o) {
        //int num1 = peptides.size();
        //int num2 = ((SimpleProteinItem)o).getPeptideItems().size();
        //double num1 = getSumZScore();
        //double num2 = ((SimpleProteinItem)o).getSumZScore();
        double num1 = getNumPeptides();
        double num2 = ((SimpleProteinItem)o).getNumPeptides();
        //double num2 = getLength();
        //double num1 = ((SimpleProteinItem)o).getLength();
        if(num1 == num2) {
            return 0;
        } else if(num1 < num2) {
            return 1;
        } else {
            return -1;
        }

    }

}

