/**
 * @file DiffMod.java
 * This is the source file for edu.scripps.pms.blindptm.Modification
 * @author Tao Xu
 * @date $Date: 2007/09/24 21:30:27 $
 */

package edu.scripps.pms.blindptm;


import java.util.*;
import java.text.DecimalFormat;


public class DiffMod {
    private double massShift;
    private char symbol;
    //private HashSet<Modification> mods = new HashSet<Modification>();
    int hashcode;
    private String info;
    private boolean [] modifiables = new boolean[256];
    private double dbinwidthMassShift;
    private double neutralLoss;
    private double dbinwidthNeutralLoss;
    private int intMassShift;    
    public static DecimalFormat threeDigits = new DecimalFormat("0.000");
//    private ArrayList<Double> differentialModifications = new ArrayList<Double>();

    public DiffMod(double massShift, char symbol) {
        this.symbol = symbol;
        this.massShift = massShift;
        intMassShift = (int)(massShift*MassCalculator.MASSACCURACYFACTOR+0.5);
        dbinwidthMassShift = massShift*MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        info = "(" + threeDigits.format(massShift) + ")";
    }
    public DiffMod(double massShift, double neutralLoss, char symbol) {
        this.symbol = symbol;
        this.massShift = massShift;
    
        dbinwidthMassShift = massShift*MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        this.neutralLoss = neutralLoss;  
        dbinwidthNeutralLoss = neutralLoss*MassSpecConstants.DBINWIDTH;
        info = "(" + threeDigits.format(massShift) + ")";
    }
    public void setModifiable(byte residue, boolean modifiable) {
        modifiables[residue] = modifiable;
    }
    public void setModifiable(boolean [] residues) {
        modifiables = residues;
    }
    public boolean isModifiable(byte residue) {
        return modifiables[residue];
    }
    
    public double getDbinwidthNeutralLoss() {
        return dbinwidthNeutralLoss;
    }
    public double getDbinwidthMassShift() {
        return dbinwidthMassShift;
    }
    public void setMassShift(double shift) {
        massShift = shift;
        intMassShift = (int)(massShift*MassCalculator.MASSACCURACYFACTOR+0.5);
        dbinwidthMassShift = massShift*MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        info = "(" + threeDigits.format(massShift) + ")";
        
    }
    public double getMassShift() {
        return massShift;
    }
    public int getIntMassShift() {
        return intMassShift;
    }
    public char getSymbol() {
        return symbol;
    }
    public String getModInfo() {
        return info;
        //return ""+symbol;
    }
    public void setModInfo(String s) {
        info = s;
    }
    public boolean equals(DiffMod o) {
        return massShift == o.massShift;
    }
    public int hashCode() {
        return hashcode; 
    }
}



