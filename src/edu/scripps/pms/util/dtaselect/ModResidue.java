package edu.scripps.pms.util.dtaselect;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class ModResidue {

    private char residue;
    private double massShift;

    public ModResidue(char residue, double massShift)
    {
        this.residue = residue;
        this.massShift = massShift;
    }

    public char getResidue() {
        return residue;
    }

    public void setResidue(char residue) {
        this.residue = residue;
    }

    public void setMassShift(double massShift) {
        this.massShift = massShift;
    }

    public double getMassShift() {
        return massShift;
    }

}
