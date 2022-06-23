/**
 * @file ModifiedPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.ModifiedPeptideHit
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.blindssptm;



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
    public void calcNumPeaksMatched () {
        getChargeNNumPeaksMatched((ppl.getChargeState()+2)/2); 
    }
    public int [] getTheorMasses() {
        //return getChargeNTheorMasses((ppl.getChargeState()+2)/2); 
        return getChargeNTheorMasses(5); 
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

        // for S-S search, consider +1 and +2 ions if false, +3 and +4 ions if true
        boolean isYModified = false; 
        boolean isBModified = false;

        for(int i = start; i < end; i++) {
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi];

            isBModified = bd == null? isBModified : true; // temp for S-S search
            isYModified = yd == null? isYModified : true;
            
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);
//System.out.println((char)seq[i] + "\t" + bMass + "\t" + yMass);
            
            int int1B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int1Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSH)/2);
            int int2Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSH)/2);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASS2H)/3);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASS2H)/3);
            int int4B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASS3H)/4);
            int int4Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASS3H)/4);
            // doubly charged fragment ions

            if(isYModified) {
                // consider +3 and +4 fragment ions                
                if(int3Y < lastTrue && int3Y > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int3Y]) numPeaksMatched++;
                }
                if(int4Y < lastTrue && int4Y > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int4Y]) numPeaksMatched++;
                }
            } else {
                // consider +1 and +2 fragment ions                

                if(int1Y < lastTrue && int1Y > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int1Y]) numPeaksMatched++;
                }
                if(int2Y < lastTrue && int2Y > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int2Y]) numPeaksMatched++;
                }
            }


            if(isBModified) {
                
                // consider +3 and +4 fragment ions                
                if(int3B < lastTrue && int3B > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int3B]) numPeaksMatched++;
                }
                if(int4B < lastTrue && int4B > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int4B]) numPeaksMatched++;
                }

            } else {
                // consider +1 and +2 fragment ions                
                if(int1B < lastTrue && int1B > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int1B]) numPeaksMatched++;
                }
                if(int2B < lastTrue && int2B > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int2B]) numPeaksMatched++;
                }

            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

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

        // for S-S search, consider +1 and +2 ions if false, +3 and +4 ions if true
        boolean isYModified = false; 
        boolean isBModified = false;

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

            isBModified = bm == null? isBModified : true; // temp for S-S search
            isYModified = ym == null? isYModified : true;
  
            if(isBModified) {

            } else {
                processSinglyChargedBIon(theorMass, bMass);
                processDoublyChargedBIon(theorMass, bMass);

            }

            if(isYModified) {

            } else {
                processSinglyChargedYIon(theorMass, yMass);

                processDoublyChargedYIon(theorMass, yMass);

            }
          

            // for highly charged fragment ions
            processHighlyChargedIons(theorMass, bMass, yMass, z);
        }
        return theorMass;
    }
    public boolean equals(ModifiedPeptideHit m) {
        return super.equals(m) && mods.equals(m.mods);
    }
}
