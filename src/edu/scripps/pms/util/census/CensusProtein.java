package edu.scripps.pms.util.census;

/**
 * @author Tao Xu    
 * @version $Id
 */
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

public class CensusProtein {

    private ArrayList<String> proteinlines = new ArrayList<String>();
    private ArrayList<String> accessions = new ArrayList<String>();
    private ArrayList<CensusPeptide> peptides = new ArrayList<CensusPeptide>(10);
    private CensusPeptide minRatioPeptide = null; 
    private CensusPeptide maxRatioPeptide = null; 

    private double ratio = 0;
    private int numPeptide = 0;
    private int spectrumCount = 0;
    private String description = null;
    private String representativeAccession = null;
    private String representativeProteinLine = null;


    public CensusProtein(String pline) {
        addProteinLine(pline); 
    }
    public void addProteinLine(String pline) {

        String [] arr = pline.split("\t");
        if(arr.length > 6) {
            ratio = Double.parseDouble(arr[2]);
            numPeptide = Integer.parseInt(arr[4]);
            spectrumCount = Integer.parseInt(arr[5]);
            description = arr[6];
            representativeAccession = arr[1].substring(0, arr[1].indexOf("."));
            representativeProteinLine = pline;
        }

        proteinlines.add(pline);

    }

    public String getDescription() {
        return description;
    }

    public String getRepresentativeAccession() {
        return representativeAccession;
    }

    public double getProteinRatio() {
        return ratio;
    }

    public int getNumPeptide() {
        return numPeptide;
    }

    public int getSpectrumCount() {
        return spectrumCount;
    }

    public void addPeptide(String sline) {
        CensusPeptide peptide = new CensusPeptide(sline);
        if(minRatioPeptide == null) {
            minRatioPeptide = peptide;
        } else {
            minRatioPeptide = minRatioPeptide.getRatio() > peptide.getRatio()? peptide : minRatioPeptide;
        }       
        
        if(maxRatioPeptide == null) {
            maxRatioPeptide = peptide;
        } else {
            maxRatioPeptide = maxRatioPeptide.getRatio() < peptide.getRatio()? peptide : maxRatioPeptide;
        }       
    }
 
    public CensusPeptide getMaxRatioPeptide() {
        return maxRatioPeptide;
    }

    public CensusPeptide getMinRatioPeptide() {
        return minRatioPeptide;
    }
}
