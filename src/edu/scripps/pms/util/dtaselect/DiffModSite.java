package edu.scripps.pms.util.dtaselect;

/**
 * @author  Tao Xu
 * @version $Id: DiffModSite.java,v 1.2 2010/10/29 04:16:16 taoxu Exp $
 */
public class DiffModSite
{
    private int site;
    private char residue;
    private double massshift;
 
    public DiffModSite(int position, double mass, char c) {
        site = position;
        massshift = mass;
        residue = c;
         
    }

    public char getResidue() {
        return residue;
    }
    public int getSite() {
        return site;
    }
    public double getMassShift() {
        return massshift;
    }
}
