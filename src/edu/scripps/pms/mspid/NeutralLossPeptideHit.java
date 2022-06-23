/**
 * @file NeutralLossPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.NeutralLossPeptideHit
 * @author Tao Xu
 * @date $Date: 2007/03/21 21:41:47 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

import java.util.*;

public class NeutralLossPeptideHit extends ModifiedPeptideHit {
     

    public NeutralLossPeptideHit(Fasta parent, int start, int end, Modifications m) {
        super(parent, start, end, m);
    }  
    public NeutralLossPeptideHit(Peptide p, Modifications m) {
        super(p, m);
    }


    public int [] getCharge2TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
  
        // here assume all the diffmods have the same neutral loss value
        int numYNeutralLoss = 0;
        int numBNeutralLoss = 0;

        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            numYNeutralLoss = processNeutralLoss(theorMass, yMass, numYNeutralLoss, ym);
            numBNeutralLoss = processNeutralLoss(theorMass, bMass, numBNeutralLoss, bm);
            
        }
        return theorMass;
    }
    private int processNeutralLoss(int [] theorMass, double mass, int nl, DiffMod m) {
        if(m == null || m.getDbinwidthNeutralLoss() == 0) {
            return nl;
        } else {
            nl++;
        }
        double neutralLoss = m.getDbinwidthNeutralLoss();
        for(int i = 1; i <= nl; i++) {
            int index = (int) (mass - i*neutralLoss);
            if(index > 0 && index < theorMass.length) {
                theorMass[index] = 50;
            }
        }
        return nl;
    }
    public int[] getCharge3TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        // here assume all the diffmods have the same neutral loss value
        int numYNeutralLoss = 0;
        int numBNeutralLoss = 0;

        for(int i = 0; i < numAA; i++) {

            int yi = yIndex -i;

            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);

            numYNeutralLoss = processNeutralLoss(theorMass, yMass, numYNeutralLoss, ym);
            numBNeutralLoss = processNeutralLoss(theorMass, bMass, numBNeutralLoss, bm);
        }
        return theorMass;
    }

}
