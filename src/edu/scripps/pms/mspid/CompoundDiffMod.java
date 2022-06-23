package edu.scripps.pms.mspid;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Titus Jung titusj@scripps.edu on 9/30/19.
 */
public class CompoundDiffMod extends DiffMod {

    protected List<Double> bIonList = null;
    protected List<Double> yIonlist = null;

    public CompoundDiffMod(double massShift, char symbol, List<Double> bIonList, List<Double> yIonlist) {
        super(massShift, symbol);
        this.bIonList = bIonList;
        this.yIonlist = yIonlist;
    }

    public List<Double> getbIonList() {
        return bIonList;
    }

    public List<Double> getyIonlist() {
        return yIonlist;
    }
}
