/**
 * @file TopdownProtein.java
 * This is the source file for edu.scripps.pms.util.spectrum.TopdownProtein
 * @author Tao Xu
 * @date $Date: 2009/08/01 00:51:55 $
 */



package edu.scripps.pms.topdown;

import java.util.ArrayList;
//import java.util.List;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.spectrum.Peak;

public class TopdownProtein extends Fasta {

    private ArrayList<Peak> theoreticalNTermPeaks;
    private ArrayList<Peak> theoreticalCTermPeaks;
    private double [] ntermTruncations;
    private double [] ctermTruncations;
    private boolean [] ntermBoolMasses;
    private boolean [] ctermBoolMasses;

    public TopdownProtein(String defline, String sequence, SearchParams sp, MassCalculator mc) {
        super(defline, sequence);
        computeTheoreticalPeaks(sp, mc);
    }
    public TopdownProtein(Fasta f, SearchParams sp, MassCalculator mc) {
        super(f.getOriginalDefline(), f.getSequence());
        computeTheoreticalPeaks(sp, mc);
    }

    private void computeTheoreticalPeaks(SearchParams sp, MassCalculator mc) {
        //fragMasses = mc.getFragMasses(); no big difference
        //double yMass = getYStart(); 
        //double bMass = getBStart();

        int numResidues = getLength();
        theoreticalNTermPeaks = new ArrayList<Peak>(numResidues);
        theoreticalCTermPeaks = new ArrayList<Peak>(numResidues);

       // maxShift = (int)(PRECISIONFACTOR*fragTolerance);

        double yMass = MassSpecConstants.MASSH2O + sp.getStaticCTermMod(); 
        double bMass = sp.getStaticNTermMod();

        int end = numResidues - 1; 
        int yIndex = end; // index for y ion, i is used for b ion

        for(int i = 0; i < end; i++) {
            bMass += mc.getFragmentMass(sequence[i]);
            theoreticalNTermPeaks.add(new Peak(bMass, 1));
            yMass += mc.getFragmentMass(sequence[yIndex--]);
            theoreticalCTermPeaks.add(new Peak(yMass, 1));
        }
        int maxntermtr = sp.getMaxNTermTruncation();
        int maxctermtr = sp.getMaxCTermTruncation();
        ntermTruncations = new double[maxntermtr];
        ctermTruncations = new double[maxctermtr];
        ntermTruncations[0] = mc.getFragmentMass(sequence[0]); 
        ctermTruncations[0] = mc.getFragmentMass(sequence[end--]); 
        for(int i = 1; i < maxntermtr && i < numResidues; i++) {
            ntermTruncations[i] = ntermTruncations[i-1] + mc.getFragmentMass(sequence[i]); 
        } 
        for(int i = 1; i < maxctermtr && i < numResidues; i++) {
            ctermTruncations[i] = ctermTruncations[i-1] + mc.getFragmentMass(sequence[end--]); 
        } 
    }
}
