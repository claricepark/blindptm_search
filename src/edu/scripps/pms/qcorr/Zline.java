/**
 * @file Zline.java
 * This is the source file for edu.scripps.pms.util.spectrum.Zline
 * @author Tao Xu
 * @date $Date: 2007/02/27 01:31:02 $
 */



package  edu.scripps.pms.qcorr;

import java.util.ArrayList;
import java.util.List;

public class Zline {
    
    private int chargeState = 1; // default of chargeState is set to 1
    private float m2z; 
    private ArrayList<String> dlines = new ArrayList<String>();
    public Zline (int chargeState, float m2z) {
        this.chargeState = chargeState;
        this.m2z = m2z;
    }

    public int getChargeState() {
        return chargeState;
    }
    public float getM2z() {
        return m2z;
    }
    public void setM2z(float m2z) {
        this.m2z = m2z;
    }
    public void addDline(String l) {
        dlines.add(l);
    }
    public List<String> getDlines() {
        return dlines;
    }
}
