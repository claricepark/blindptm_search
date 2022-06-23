/**
 * @file MassSpecConstants.java
 * This is the source file for edu.scripps.pms.util.spectrum.MassSpecConstants
 * @author Tao Xu
 * @date $Date
 */

package edu.scripps.pms.blindssptm;


public class MassSpecConstants {

    public static float MASSH = 1.007825f;
    public static float MASS2H = 2*1.007825f;
    public static float MASS3H = 3*1.007825f;
    public static float MASS4H = 4*1.007825f;
    public static final float log10 = 0.434294481f;
    public static final float MASSH2O = 15.9949145f + 2*MASSH;
    public static final float MASSH2 = MASSH * 2;
    public static final float MASSH3O = MASSH2O + MASSH;
    public static final float DBINWIDTH = 1.f/1.0005079f;
    public static final float MASSHDB = MASSH*DBINWIDTH;
    public static final float MASSH3ODB = MASSH3O*DBINWIDTH;
    public static final String AVGISOTOPE = "avg";
    public static final String MONOISOTOPE = "mono";
    // Thermo source is about 1.002806, which does not look right
    public static final float MASSDIFFC12C13 = 1.003354826f;
    public static final float MASSDIFFN14N15 = 0.997034968f;
}
