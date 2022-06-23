/**
 * @file TerminalModification.java
 * This is the source file for edu.scripps.pms.blindptm.Modification
 * @author Tao Xu
 * @date $Date: 2007/06/07 22:43:48 $
 */

package edu.scripps.pms.blindptm;


import java.util.*;


public class TerminalModification {
    private double massShift;
    private char symbol;
    private String info;
 
    public TerminalModification(char symbol, double massShift) {
        this.symbol = symbol;
        this.massShift = massShift;
        info = "(" + massShift + ")";
    }
    public double getMassShift() {
        return massShift;
    }
    public char getSymbol() {
        return symbol;
    }
    public String getModInfo() {
        return info;
    }
    public void setModInfo(String s) {
        info = s;
    }
    public String getInfoAsSymbol() {
        return "" + symbol;
    }

    public String getInfoAsMassShift() {
        return "(" + massShift + ")";
    }
}



