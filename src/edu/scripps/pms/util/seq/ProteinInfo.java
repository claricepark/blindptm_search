package edu.scripps.pms.util.seq;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.FastaReader;
//import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.dtaselect.Protein;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// 
public class ProteinInfo {


    /* Loop through each residue in the database sequence for this
     * locus and accumulate molecular weight.  Use a simple algorithm
     * to estimate pI.
     */
    public static float [] calculatePI(String Sequence) {
	//String         Sequence = this.TrimmedSequence();
	int            length = Sequence.length();
	int            looper;
	int            countLys = 0;
	int            countArg = 0;
	int            countHis = 0;
	int            countAsp = 0;
	int            countGlu = 0;
	int            countCys = 0;
	int            countTyr = 0;
	float          currentPH = 7.0f;
	float          currentJump = 3.5f;
	float          currentCharge;
	float          lastCharge = 0;
	float          mwaccum = 0f;
	char           currentResidue;
	if (length > 0) {
	    for (looper = 0; looper < length; looper++) {
		currentResidue = Sequence.charAt(looper);
		switch (currentResidue) {
		case 'A':
		    mwaccum += 71.0;
		    break;
		case 'C':
		    mwaccum += 103.0;
		    countCys ++;
		    break;
		case 'D':
		    mwaccum += 115.0;
		    countAsp ++;
		    break;
		case 'E':
		    mwaccum += 129.0;
		    countGlu ++;
		    break;
		case 'F':
		    mwaccum += 147.0;
		    break;
		case 'G':
		    mwaccum += 57.0;
		    break;
		case 'H':
		    mwaccum += 137.0;
		    countHis ++;
		    break;
		case 'I':
		    mwaccum += 113.0;
		    break;
		case 'K':
		    mwaccum += 128.0;
		    countLys ++;
		    break;
		case 'L':
		    mwaccum += 113.0;
		    break;
		case 'M':
		    mwaccum += 131.0;
		    break;
		case 'N':
		    mwaccum += 114.0;
		    break;
		case 'P':
		    mwaccum += 97.0;
		    break;
		case 'Q':
		    mwaccum += 128.0;
		    break;
		case 'R':
		    mwaccum += 156.0;
		    countArg ++;
		    break;
		case 'S':
		    mwaccum += 87.0;
		    break;
		case 'T':
		    mwaccum += 101.0;
		    break;
		case 'V':
		    mwaccum += 99.0;
		    break;
		case 'W':
		    mwaccum += 186.0;
		    break;
		case 'Y':
		    mwaccum += 176.0;
		    countTyr ++;
		    break;
		}
	    }
	    /* Use a bracketing strategy to identify the isoelectric
	     * point.  Calculate charge at pH of 7, and then move up
	     * 3.5 or down 3.5 as necessary.  Make each successive
	     * move up or down only half as large.  Keep going until
	     * two successive charges reported match to one place past
	     * the decimal point.
	     */
	    currentCharge = chargeAtPh(currentPH, countLys,
					       countArg, countHis,
					       countAsp, countGlu,
					       countCys, countTyr);
	    while (roundTo(currentCharge,1) != roundTo(lastCharge,1)) {
		//		System.out.println("pH:\t" + new Float(currentPH).toString() 
		//              + "\tCharge\t" + new Float(currentCharge).toString());
		if (currentCharge > 0)
		    currentPH += currentJump;
		else
		    currentPH -= currentJump;
		currentJump /= 2;
		lastCharge = currentCharge;
		currentCharge = chargeAtPh(currentPH,
						   countLys, countArg,
						   countHis, countAsp,
						   countGlu, countCys,
						   countTyr);
		if ( (currentPH > 14) || (currentPH < 0) ) {
		    System.out.println("pI can't be figured for " + Sequence);
		    System.exit(0);
		}
	    }
	}
        float [] result = new float[2];
        result[0] = currentPH;
        result[1] = mwaccum;
        return result;
	//return currentPH;

    }

    // Because I don't like the Math.roundTo function, I wrote this one.
    public static float roundTo(float Value, int Places) {
	//Converts a value to a rounded value
	if (Places == 0) {
	    return new Integer(Math.round(Value)).intValue();
	}
	else {
	    double Multiplier = Math.pow(10,Places);
	    return new Double(Math.rint(Value*Multiplier)/Multiplier).floatValue();
	}
    }
    
    /* Determine the charge on a theoretical protein at a given pH
       given its relevant composition */
    public static float chargeAtPh(float pH, int countLys, int countArg, int countHis,
				   int countAsp, int countGlu, int countCys, int countTyr) {
	//Start out accumulator with charge of termini
	float          accum = percentPositive(pH, 8.0f) - percentNegative(pH, 3.1f);
	accum += countLys * percentPositive(pH, 10.0f);
	accum += countArg * percentPositive(pH, 12.0f);
	accum += countHis * percentPositive(pH, 6.5f);
	accum -= countAsp * percentNegative(pH, 4.4f);
	accum -= countGlu * percentNegative(pH, 4.4f);
	accum -= countCys * percentNegative(pH, 8.5f);
	accum -= countTyr * percentNegative(pH, 10.0f);
	return accum;
    }

    /* What percentage of ions at given pK at given pH would be
       negatively charged? */
    public static float percentNegative(float pH, float pK) {
	double         concentrationRatio = Math.pow(10,pH - pK);
	return new Double(concentrationRatio / (concentrationRatio + 1)).floatValue();
    }
    /* What percentage of ions of given pK at given pH would be
       postively charged? */
    public static float percentPositive(float pH, float pK) {
	double         concentrationRatio = Math.pow(10f,pK - pH);
	return new Double(concentrationRatio / (concentrationRatio + 1)).floatValue();
    }
    public static float calculateKyteDoolittle(String Sequence) {
        int            length = Sequence.length();
        int            looper;
        float          mwaccum = 0f;
        char           currentResidue;
        if (length > 0) {
            for (looper = 0; looper < length; looper++) {
                currentResidue = Sequence.charAt(looper);
                switch (currentResidue) {
                case 'A':
                    mwaccum += 18.0;
                    break;
                case 'C':
                    mwaccum += 25.0;
                    break;
                case 'D':
                    mwaccum -= 35.0;
                    break;
                case 'E':
                    mwaccum -= 35.0;
                    break;
                case 'F':
                    mwaccum += 28.0;
                    break;
                case 'G':
                    mwaccum -= 4.0;
                    break;
                case 'H':
                    mwaccum -= 32.0;
                    break;
                case 'I':
                    mwaccum += 45.0;
                    break;
                case 'K':
                    mwaccum -= 39.0;
                    break;
                case 'L':
                    mwaccum += 38.0;
                    break;
                case 'M':
                    mwaccum += 19.0;
                    break;
                case 'N':
                    mwaccum -= 35.0;
                    break;
                case 'P':
                    mwaccum -= 16.0;
                    break;
                case 'Q':
                    mwaccum -= 35.0;
                    break;
                case 'R':
                    mwaccum -= 45.0;
                    break;
                case 'S':
                    mwaccum -= 8.0;
                    break;
                case 'T':
                    mwaccum -= 7.0;
                    break;
                case 'V':
                    mwaccum += 42.0;
                    break;
                case 'W':
                    mwaccum -= 9.0;
                    break;
                case 'Y':
                    mwaccum -= 13.0;
                    break;
                }
            }
        }
        return mwaccum/10.0f;
    }

    public static int getNumTrypticPeptides(Fasta f, Protease p, MassCalculator mc, double highLimit, double lowLimit, int minlength) {
        double [] precMasses = mc.getPrecMasses();
        int counts = 0;
        int proteinlength =  f.getLength();
        int lastIndex = f.getLength() - 1;
        byte [] seq = f.getSequenceAsBytes();
            // mass of H2O need to be added
            // because the mass of the residues does not count H2O
        double mass = MassSpecConstants.MASSH3O;
        double tempmass = 0;
        int templeft = 0;
        //int rightEnd = -1;
        int leftEnd = 0;

        while(leftEnd < proteinlength) {

            int rightEnd = leftEnd;
            
            while(rightEnd < lastIndex) {
                rightEnd++;
                if(rightEnd == lastIndex || seq[rightEnd-1] == 'K' ||  seq[rightEnd-1] == 'R') { // tryptic end
                    if(rightEnd - leftEnd > minlength-1) {
                        double precmass = mc.getPrecursorMass(seq, leftEnd, rightEnd);
                        if(precmass >= lowLimit && precmass <= highLimit) {

                            counts++;
                        }
                    }
                    break;
                }
            }
            leftEnd = rightEnd+1;

        }

        return counts;

    }

}

