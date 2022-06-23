/**
 * @file MassSpecConstants.java
 * This is the source file for edu.scripps.pms.util.spectrum.MassSpecConstants
 * @author Tao Xu
 * @date $Date: 2009/02/02 05:11:43 $
 */

package edu.scripps.pms.mspid;

import java.text.DecimalFormat;
public class MassSpecConstants {

    public static float MASSH = 1.007825f;
    public static float MASSPROTON = 1.007276466f;

    public static final float DBINWIDTH = 1.f/1.0005079f;

    public static final float log10 = 0.434294481f;
    public static final float MASSO = 15.9949145f;
    public static final float MASSC = 12.0f;
    public static final float MASSCO = (MASSC+MASSO);

    public static final float MASSCOD2 = MASSCO/2f;

    public static final float MASSN = 14.003074008f;
    public static final float MASSNH3 = (3*MASSH+MASSN);
    public static final float MASSNH3D2  = MASSNH3/2f;

    public static final float MASSH2O = (15.9949145f + 2*MASSH);
    public static final float MASSH20D2 = MASSH2O/2.f;
    //public static final float MASSH2 = MASSH * 2;
    // MASSH2O + MASSPROTON
    public static final float MASSH3O = MASSH2O + MASSPROTON;
    public static final float MASSHDB = MASSH*DBINWIDTH;
    public static final float MASSPROTONDB = MASSPROTON*DBINWIDTH;
    public static final float MASSH3ODB = MASSH3O*DBINWIDTH;
    public static final String AVGISOTOPE = "avg";
    public static final String MONOISOTOPE = "mono";
    // Thermo source is about 1.002806, which does not look right
    public static final float MASSDIFFC12C13 = 1.003354826f;
    public static final float MASSDIFFC12C13DB = 1.003354826f*DBINWIDTH;
    public static final float MASSDIFFN14N15 = 0.997034968f;
    public static final int NUMCHARS = 256;

    public static final DecimalFormat threeDigits = new DecimalFormat("0.000");
    public static final DecimalFormat fourDigits = new DecimalFormat("0.0000");
    public static final DecimalFormat fiveDigits = new DecimalFormat("0.00000");
    public static final DecimalFormat twoDigits = new DecimalFormat("0.00");
    public static final DecimalFormat oneDigit = new DecimalFormat("0.0");
}
