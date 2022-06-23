/**
 * @file DiffModSite.java
 * This is the source file for edu.scripps.pms.mspid.Modification
 * @author Tao Xu
 * @date $Date: 2007/04/12 22:49:56 $
 */

package edu.scripps.pms.topdown;


import java.util.*;


public class DiffModSite {
//    private ArrayList<Double> differentialModifications = new ArrayList<Double>();
    private int position;
    private HashSet<DiffMod> mods;

    public DiffModSite(int pos, HashSet<DiffMod> mods) {
        position = pos;
        this.mods = mods;
    }
    public int numDiffMods() {
        return mods.size();
    }
    public Iterator<DiffMod> getDiffMods() {
        return mods.iterator();
    }
    public int getPosition() {
        return position;
    }
}



