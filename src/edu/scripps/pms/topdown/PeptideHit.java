/**
 * @file PeptideHit.java
 * This is the source file for edu.scripps.pms.topdown.PeptideHit
 * @author Tao Xu
 * @date $Date: 2007/12/02 18:42:03 $
 */



package edu.scripps.pms.topdown;

import edu.scripps.pms.util.seq.Fasta;

public class PeptideHit implements Comparable<PeptideHit> {

    // description line of this PeptideHit sequence
    // the sequence string of this PeptideHit
    
    //private int chargeState;
    private int numPeaksMatched = 0;
    private int numTheroticPeaks = 0;
    private double probability;
    protected static double [] masses;
    protected ProcessedPeakList ppl;    
    protected Peptide peptide;
    public static final boolean isModified = false;
    public static final int maxNumShift = 1000; // max addition or deletion or mass 1000 from both nterm and cterm
    public static final int maxNumShift2 = maxNumShift*2;

    private int [] numNTermPeaksMatched = new int[maxNumShift2]; // n term
    private int [] numCTermPeaksMatched = new int[maxNumShift2]; // c term
    
    private int ntermIndex; // b ions
    private int ctermIndex; // y ions 
    private double ntermMassShift; // b ions
    private double ctermMassShift; // y ions

    //private TopdownProtein;
    private int [][] ntermMatchedPeaks; // first dimension for truncation and the second dimension for modification
    private int [][] ctermMatchedPeaks;
    private int [] ntermTestedPeaks; // one element for each 
    private int [] ctermTestedPeaks; // 

    public double getNTermMassShift() {
        return ntermMassShift;
    }
    public double getCTermMassShift() {
        return ctermMassShift;
    }
    public PeptideHit(Fasta parent, int start, int end) {
        peptide = new Peptide(parent, start, end); 
    }  
    public PeptideHit(Peptide p) {
        peptide = p;
    }
    public void setProcessedPeakList(ProcessedPeakList ppl) {
        this.ppl = ppl;
        masses = ppl.getFragMasses();
//for(double f : masses) {
//    if(f != 0)
//    System.out.println(f + "\t");
//}
    }
    public boolean isModified() {
        return isModified;
    }
    
    public void calcNumPeaksMatched () {
        double maxmz = ppl.getPeakList().getMaxM2z();
        double minmz = ppl.getPeakList().getMinM2z();
//System.out.println("maxmz: " + maxmz + "\tminmz: " + minmz);
        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getYStart(); // from super class
        double bMass = getBStart();  // from super class
        
        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion
        int maxIndex = boolMasses.length - maxNumShift;
         
        for(int i = start; i < end; i++) {
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            
            //if(bMass >= minmz && bMass <= maxmz) {
            if(tempB >= maxNumShift && bMass <= maxmz && tempB < maxIndex) {
                numTheroticPeaks++;
                checkMatches(boolMasses, tempB, numNTermPeaksMatched); 
            }
            //if(yMass >= minmz && yMass <= maxmz) {
            if(tempY > maxNumShift && yMass < maxmz && tempY < maxIndex) {
                numTheroticPeaks++;
                checkMatches(boolMasses, tempY, numCTermPeaksMatched); 
            }
            if(bMass >= maxmz && yMass >= maxmz) {
                break;
            }
        }
        ntermIndex = getMaxShiftIndex(numNTermPeaksMatched); 
        ctermIndex = getMaxShiftIndex(numCTermPeaksMatched);
        ntermMassShift = (ntermIndex-maxNumShift)/(ppl.PRECISIONFACTOR+0.0); 
        ctermMassShift = (ctermIndex-maxNumShift)/(ppl.PRECISIONFACTOR+0.0); 
        setNumPeaksMatched(ppl.getPTrue());
    }
    private int getMaxShiftIndex(int [] numPeaksMatched) {
        int maxindex = 0;
        int maxvalue = 0;
        for(int i = 0; i < maxNumShift2; i++) {
            if(numPeaksMatched[i] > maxvalue) {
                maxindex = i;
                maxvalue = numPeaksMatched[i];
            }
        }
//System.out.println("maxindex: " + maxindex + "\tnumPeaksMatched: " + numPeaksMatched[maxindex]); 
        return maxindex;
    }
    private void checkMatches(boolean [] boolMasses, int mass, int [] numPeaksMatched) {
         
        for(int i = 0; i < maxNumShift2; i++) {
            int index = mass + i - maxNumShift;
              
            //if(index > 0 && index < boolMasses.length && boolMasses[index]) {
            if(boolMasses[index]) {
                numPeaksMatched[i]++;
            }
        } 
    }
    public double getTheorMass() {
        double mass = ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +
                     ppl.getSearchParams().getStaticTerminalMods();
   
        return mass;
    }
    
    public double getDbinwidthBStart() {
        return ppl.getDbinwidthBStart();
    }
    public double getDbinwidthYStart() {
        return ppl.getDbinwidthYStart();
    }
    public double getBStart() {
        return ppl.getBStart();
    }
    public double getYStart() {
        return ppl.getYStart();
    }

    public int [] getTheorMasses() {
        double maxmz = ppl.getPeakList().getMaxM2z();
        double minmz = ppl.getPeakList().getMinM2z();
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        //double bMass = MassSpecConstants.MASSHDB + 0.5f;
        //double yMass = MassSpecConstants.MASSH3ODB + 0.5f;
        double bMass = getDbinwidthBStart() + 0.5f;
        double yMass = getDbinwidthYStart() + 0.5f;
        //double tempb = 0;
        //double tempy = 0;
//if(bMass > 20)
//System.out.println("ystart: " + yMass + "\tbstart: " + bMass);

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            bMass += masses[seq[i]];
            yMass += masses[seq[yIndex-i]];
            if(bMass >= minmz && bMass <= maxmz) {
                processBIon(theorMass, bMass);
            }

            if(yMass >= minmz && yMass <= maxmz) {
                processYIon(theorMass, yMass);
            }
            if(bMass >= maxmz && yMass >= maxmz) {
                break;
            }
        }
        return theorMass;
    }

    protected void processBIon(int [] theorMass, double mass) {
//System.out.println("in singlyChargedBIon, mass: " + mass);
        try {
            int intMass = (int)mass;
            theorMass[intMass] = 50;
            int index = intMass + 1;
            if(theorMass[index] < 25) theorMass[index] = 25;

            index++; // -= 2; 
            if(theorMass[index] < 25) theorMass[index] = 25;

            index = intMass - 18; // H2O loss
            if(theorMass[index] < 10) theorMass[index] = 10;
            index = intMass - 28; // CO loss
            if(theorMass[index] < 10) theorMass[index] = 10;
            index = intMass - 17; // NH3 loss
            if (theorMass[index] < 10) theorMass[index] = 10;
        } catch(Exception e) {}// igore exception caused by weird aa residue
       
    }
    protected void processYIon(int [] theorMass, double mass) {
//System.out.println("in singlyChargedYIon, mass: " + mass);
        try { 
            int intMass = (int)mass;
            //System.out.println("intmass: " + intMass);
            theorMass[intMass] = 50;
            int index = intMass + 1;
            if(theorMass[index] < 25)  theorMass[index] = 25;
 
            index++; // -= 2; 
            if(theorMass[index] < 25) theorMass[index] = 25;

            index = intMass - 17; // loss NH3
            if(theorMass[index] < 10) theorMass[index] = 10;
        
        } catch(Exception e) {} 
    }
    public boolean equals(PeptideHit p) {
        
        if(p != null) {
            return getSequence().equals(p.getSequence()); 
        }
        return false;
    }
    public int compareTo(PeptideHit p) {

        return p.numPeaksMatched - this.numPeaksMatched; 
        /*
        if(p.probability < this.probability) {
            return 1;
        } else if(p.probability > this.probability) {
            return -1;
        } else { return 0; }
       */
    }
    public void setNumPeaksMatched(double ptrue) {
        numPeaksMatched = numCTermPeaksMatched[ctermIndex] + numNTermPeaksMatched[ntermIndex];
        probability = 1.0/numPeaksMatched;
        //probability = DistributionCalculator.getBinomialSum(ptrue, numTheroticPeaks, numPeaksMatched);
//System.out.println("Probability: " + probability + "\t" + numTheroticPeaks + "\t" + numPeaksMatched );
    }
    public Fasta getParent() {
        return peptide.getParent();
    }
    public double getProbability() {
        return probability;
    }
    public String getExtendedSequence() {  
        return peptide.getExtendedSequence();
    }
    public String getSequence() {
        return peptide.getSequence(); 
    }
    public int getStart() {
        return peptide.getStart();
    }
    public String getDefline() {
        return peptide.getDefline();
    }
    public int getEnd() {
        return peptide.getEnd();
    }

    public int getNumPeaksMatched() {
        return numPeaksMatched;
    }
    public int getLength() {
        return peptide.getLength();
    }
    public String getAccession() {
        return peptide.getAccession();
    }
    public String getDescription() {
        return peptide.getDescription();
    }
    public int getNumPeaks() {
        //return 2*2*(chargeState-1)*(end - start); // 2*(lengh - 1)
        return numTheroticPeaks;
    }
}
