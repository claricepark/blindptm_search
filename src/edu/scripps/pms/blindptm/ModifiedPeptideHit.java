/**
 * @file ModifiedPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.ModifiedPeptideHit
 * @author Tao Xu
 * @date $Date: 2007/09/24 21:30:27 $
 */



package edu.scripps.pms.blindptm;

import edu.scripps.pms.util.seq.Fasta;

import java.util.*;

public class ModifiedPeptideHit extends PeptideHit {
     
    protected Modifications mods;
    protected String extendedSeq = null;
    protected String exactSeq = null; // same as extendedSeq except leading and tailing residue
    protected DiffMod [] diffMods;
    protected ArrayList<ModifiedPeptideHit> candidates;
    protected ArrayList<DiffModSite> sites = new ArrayList<DiffModSite>();
    protected int numPotentialSites = -1;
    protected double internalDiffModMassAdded = 0;
    protected int numInternalDiffModsAdded = 0;
    public static final boolean isModified = true;


    public ModifiedPeptideHit(Fasta parent, int start, int end, Modifications m) {
        super(parent, start, end);
        this.mods = m;
        diffMods = new DiffMod[end-start+1];
    }  
    public ModifiedPeptideHit(Peptide p, Modifications m) {
        super(p);
        mods = m;
        diffMods = new DiffMod[p.getEnd() - p.getStart() + 1];

    }
    public ModifiedPeptideHit(Peptide p) {
        super(p);
        mods = new Modifications();
        diffMods = new DiffMod[p.getEnd() - p.getStart() + 1];

    }
    public boolean isModified() {
        return isModified;
    }
    public void setDiffMod(int index, DiffMod diffmod) {
        diffMods[index] = diffmod;
        if(diffmod != null) {
            internalDiffModMassAdded += diffmod.getMassShift();
            numInternalDiffModsAdded++;
        }
    }
    public boolean isGoodCandidate() {
        return internalDiffModMassAdded == mods.getInternalDiffModsMass() &&
              numInternalDiffModsAdded == mods.getNumDiffMods();
    }
    public void setDiffMods(DiffMod [] dms) {
        this.diffMods = dms;
    }
    protected void generateCandidates() {
        candidates = new ArrayList<ModifiedPeptideHit>();
        calcModificationSites();
//System.out.println("numSites: " + numPotentialSites + "\t numDiffMods: " + mods.getNumDiffMods());
        if(getNumPotentialSites() < mods.getNumDiffMods()) {
            return;
        }
        candidates.add(this.copy()); 
        int numSites = sites.size(); 

//System.out.print("numPotentialSitesInPeptide: " + numPotentialSites +  "\t numDiffMods: " + mods.getNumDiffMods());
        // may have more efficient way if numSites == mods.getNumDiffMods()
        for(int i = 0; i < numSites; i++) {
           
            DiffModSite dms = sites.get(i); 
            processDiffModSite(dms);
        }
        getFinalCandidates();
    }
    private void getFinalCandidates() {
        
        ArrayList<ModifiedPeptideHit> temp = candidates;
        candidates = new ArrayList<ModifiedPeptideHit>();
        for(Iterator<ModifiedPeptideHit> it = temp.iterator(); it.hasNext();) {
            ModifiedPeptideHit m = it.next();
            if(m.isGoodCandidate()) {
                candidates.add(m);
            }
        } 
//System.out.println("\tNumTotal Candidates: " + temp.size() + "\tNumGoodCandidates: " + candidates.size());
    }
    private void processDiffModSite(DiffModSite dms) {
//System.out.println("numDiffModsForThisSite: " + dms.numDiffMods());
        int num = candidates.size();
        int position = dms.getPosition();
        for(int i = 0; i < num; i++) {
            ModifiedPeptideHit m = candidates.get(i);
            for(Iterator<DiffMod> it = dms.getDiffMods(); it.hasNext();) {
                DiffMod d = it.next();
                if(m.numInternalDiffModsAdded+1 <= mods.getNumDiffMods() &&
                          m.internalDiffModMassAdded+d.getMassShift() <= mods.getInternalDiffModsMass()) {
                    ModifiedPeptideHit newm = m.copy();
                    newm.setDiffMod(position, d);
                    candidates.add(newm); 
                }
            }
        }
    }   
    private ModifiedPeptideHit copy() {
        ModifiedPeptideHit newm = new ModifiedPeptideHit(peptide, mods);
        for(int i = 0; i < newm.diffMods.length; i++) {
            newm.setDiffMod(i, diffMods[i]);
        }              
        //newm.internalDiffModMassAdded = internalDiffModMassAdded;
        return newm;
    }
    protected void calcModificationSites() {
        if(numPotentialSites > -1) {
            return;
        }
        numPotentialSites = 0;
        int start = getStart();
        int end = getEnd();
        byte [] bytes = getParent().getSequenceAsBytes();
        for(int i = start; i <= end; i++) {
            if(mods.isModifiable(bytes[i])) {
                numPotentialSites++;
                sites.add(new DiffModSite(i - start, mods.residue2DiffMods(bytes[i])));
            }
        }
    }
    public int getNumPotentialSites() {
        if(numPotentialSites > -1) {
            calcModificationSites();
        }
        return numPotentialSites; 
    }    

    public Iterator<ModifiedPeptideHit> getAllModifiedPeptideHits() {
        generateCandidates(); 
        return candidates.iterator();
    }
    public double getTheorMass() {
       
        double m = ppl.getMassCalculator().getPrecursorMass(peptide.getSequence());
        m += internalDiffModMassAdded == 0? mods.getMassShift() : internalDiffModMassAdded;
          
        //return ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +  mods.getMassShift();
        return m;
    }
    // override the function in super class
    public double getNTermStart() { // b == n term
        if(mods != null) {
            return ppl.getNTermStart() + mods.getDiffNTermMod();
        } else {
            return ppl.getNTermStart();

        }
    }    
    // override the function in super class
    public double getNTermDbinwidthStart() { // b == n term
        if(mods != null) {
            return (ppl.getNTermStart() + mods.getDiffNTermMod())*MassSpecConstants.DBINWIDTH;
        } else {
            return ppl.getNTermStart()*MassSpecConstants.DBINWIDTH;
        }
    }    
    public double getCTermStart() { // y == c term
//System.out.println("in getCTermStart in ModifiedPeptideHit, " + mods.getDiffCTermMod());
        if(mods != null) {
            return ppl.getCTermStart() + mods.getDiffCTermMod();
        } else {     
            return ppl.getCTermStart();
        }
    }    
    public double getCTermDbinwidthStart() { // y == c term
        if(mods != null) {
            return (ppl.getCTermStart() + mods.getDiffCTermMod())*MassSpecConstants.DBINWIDTH;
        } else {

            return ppl.getCTermStart()*MassSpecConstants.DBINWIDTH;
        }
    }    

    // return sequence without leading and tailing residues
    public String getExactSequence() {
        if(exactSeq != null) {
            return exactSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        int start = peptide.getStart();
        int end = peptide.getEnd();

        double nTermMod = mods == null? 0 : mods.getDiffNTermMod();
        double cTermMod = mods == null? 0 : mods.getDiffCTermMod();
        if(nTermMod != 0) {
            sb.append("(" + nTermMod + ")"); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != 0) {
            sb.append("(" + cTermMod + ")"); 
        }

        exactSeq = sb.toString(); 
        
        return exactSeq;
    }

    public String getExtendedSequence() {
        if(extendedSeq != null) {
            return extendedSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        // 4 more bytes for the heading and tailing aa and '.'
        int len = peptide.getLength() + 4;
        int start = peptide.getStart();
        int end = peptide.getEnd();
        if (start == 0) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(start-1));
        }
        sb.append('.');

        TerminalModification nTermMod = mods.getNTerm();
        TerminalModification cTermMod = mods.getCTerm();
        if(nTermMod != null) {
            //sb.append("(" + nTermMod + ")"); 
            sb.append( nTermMod.getModInfo()); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != null) {
            //sb.append("(" + cTermMod + ")"); 
            sb.append(cTermMod.getModInfo()); 
        }
        sb.append('.');

        int lastIndex = peptide.getParent().getLength() - 1;
        if (end == lastIndex) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(end+1));
        }

        extendedSeq = sb.toString(); 
        
        return extendedSeq;
    }


    public void getChargeNNumPeaksMatched(int z) {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function 
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;
        for(int i = start; i < end; i++) {


            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);
//System.out.println((char)seq[i] + "\t" + bMass + "\t" + yMass);
            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSH)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSH)/2f);
            // doubly charged fragment ions
            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;
 
            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                numTheroticPeaks += 2;
                if(boolMasses[int2B]) numPeaksMatched++;
                if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                numTheroticPeaks += 2;
                if(boolMasses[int2Y]) numPeaksMatched++;

                if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
           // if(boolMasses[int2B++]) numPeaksMatched++;
           // if(boolMasses[int2Y++]) numPeaksMatched++;
            //numTheroticPeaks += 6;
            for(int j = 3; j < z; j++) {
                int nb = (int)(0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSH)/j);
                int ny = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSH)/j);
                if(nb < lastTrue && nb > firstTrue) { 
                    numTheroticPeaks++;
                    if(boolMasses[nb]) numPeaksMatched++;
                }
                if(ny < lastTrue && ny > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[ny]) numPeaksMatched++;
                }
            }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }
    public void getCharge3NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion
        

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSH)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSH)/2f);

            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;

            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                numTheroticPeaks += 2;
                if(boolMasses[int2B]) numPeaksMatched++;
                if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                numTheroticPeaks += 2;
                if(boolMasses[int2Y]) numPeaksMatched++;

                if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }

    public void getCharge2NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getCTermStart(); // from super class
        double bMass = getNTermStart();  // from super class

        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);

            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(tempB < lastTrue && tempB > firstTrue) {
                numTheroticPeaks += 2;
                if(boolMasses[tempB]) numPeaksMatched++;
                if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                numTheroticPeaks += 2;
                if(boolMasses[tempY]) numPeaksMatched++;

                if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }

    public int [] getCharge2TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
       
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);
        }
        return theorMass;
    }
    public int[] getCharge3TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex -i;

            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];
     
            
            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);
        }
        return theorMass;
    }

    public int [] getChargeNTheorMasses(int z) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;


        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex-i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();


            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);

            // for highly charged fragment ions
            processHighlyChargedIons(theorMass, bMass, yMass, z);
        }
        return theorMass;
    }
    public boolean equals(ModifiedPeptideHit m) {
        return super.equals(m) && mods.equals(m.mods);
    }
}
