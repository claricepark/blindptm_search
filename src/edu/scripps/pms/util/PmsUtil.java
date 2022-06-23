package edu.scripps.pms.util;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version $Id
 */


public class PmsUtil {

    /* Loop through each residue in the database sequence for this
     * locus and accumulate molecular weight.  Use a simple algorithm
     * to estimate pI.
     */
    public static float calcPi(String Sequence) {
	// Get the default amino acid masses...
	int            Length = Sequence.length();
	int            Looper;
	int            CountLys = 0;
	int            CountArg = 0;
	int            CountHis = 0;
	int            CountAsp = 0;
	int            CountGlu = 0;
	int            CountCys = 0;
	int            CountTyr = 0;
	float          CurrentPH = 7.0f;
	float          CurrentJump = 3.5f;
	float          CurrentCharge;
	float          LastCharge = 0;
	float          MWAccum = 18.0f;
	char           CurrentResidue;
	if (Length > 0) {
	    for (Looper = 0; Looper < Length; Looper++) {
		CurrentResidue = Character.toUpperCase(Sequence.charAt(Looper));
	//	MWAccum += SequestParams.AvgMasses[CurrentResidue - 65];
		switch (CurrentResidue) {
		case 'C':
		    CountCys ++;
		    break;
		case 'D':
		    CountAsp ++;
		    break;
		case 'E':
		    CountGlu ++;
		    break;
		case 'H':
		    CountHis ++;
		    break;
		case 'K':
		    CountLys ++;
		    break;
		case 'R':
		    CountArg ++;
		    break;
		case 'Y':
		    CountTyr ++;
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
	    CurrentCharge = ChargeAtPH(CurrentPH, CountLys,
					       CountArg, CountHis,
					       CountAsp, CountGlu,
					       CountCys, CountTyr);
	    while (RoundTo(CurrentCharge,1) != RoundTo(LastCharge,1)) {
		//		System.out.println("pH:\t" + new Float(CurrentPH).toString() 
		//              + "\tCharge\t" + new Float(CurrentCharge).toString());
		if (CurrentCharge > 0)
		    CurrentPH += CurrentJump;
		else
		    CurrentPH -= CurrentJump;
		CurrentJump /= 2;
		LastCharge = CurrentCharge;
		CurrentCharge = ChargeAtPH(CurrentPH,
						   CountLys, CountArg,
						   CountHis, CountAsp,
						   CountGlu, CountCys,
						   CountTyr);
		if ( (CurrentPH > 14) || (CurrentPH < 0) ) {
		    System.out.println("!!! pI can't be figured out !!!");
		    System.exit(0);
		}
	    }
	    // pI = CurrentPH; // pI
//	    MolWt = MWAccum;
	}
        return CurrentPH;
    }
    // Because I don't like the Math.roundTo function, I wrote this one.
    public static float RoundTo(float Value, int Places) {
	//Converts a value to a rounded value
	if (Places == 0) {
	    return new Integer(Math.round(Value)).intValue();
	}
	else {
	    double Multiplier = Math.pow(10,Places);
	    return new Double(Math.rint(Value*Multiplier)/Multiplier).floatValue();
	}
    }
    /* What percentage of ions of given pK at given pH would be
       postively charged? */
    public static float PercentPositive(float pH, float pK) {
	double         ConcentrationRatio = Math.pow(10f,pK - pH);
	return new Double(ConcentrationRatio / (ConcentrationRatio + 1)).floatValue();
    }
    /* What percentage of ions at given pK at given pH would be
       negatively charged? */
    public static float PercentNegative(float pH, float pK) {
	double         ConcentrationRatio = Math.pow(10,pH - pK);
	return new Double(ConcentrationRatio / (ConcentrationRatio + 1)).floatValue();
    }
    /* Determine the charge on a theoretical protein at a given pH
       given its relevant composition */
    public static float ChargeAtPH(float pH, int CountLys, int CountArg, int CountHis,
				   int CountAsp, int CountGlu, int CountCys, int CountTyr) {
	//Start out accumulator with charge of termini
	float          Accum = PercentPositive(pH, 8.0f) - PercentNegative(pH, 3.1f);
	Accum += CountLys * PercentPositive(pH, 10.0f);
	Accum += CountArg * PercentPositive(pH, 12.0f);
	Accum += CountHis * PercentPositive(pH, 6.5f);
	Accum -= CountAsp * PercentNegative(pH, 4.4f);
	Accum -= CountGlu * PercentNegative(pH, 4.4f);
	Accum -= CountCys * PercentNegative(pH, 8.5f);
	Accum -= CountTyr * PercentNegative(pH, 10.0f);
	return Accum;
    }
}
