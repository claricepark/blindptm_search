/**
 * @file TerminalModification.java
 * This is the source file for edu.scripps.pms.topdown.Modification
 * @author Tao Xu
 * @date $Date: 2007/04/12 22:49:56 $
 */

package edu.scripps.pms.topdown;


import java.util.*;


public class TerminalModification {
    private double massShift;
    private char symbol;
    private boolean isStatic;
 
    public TerminalModification(char symbol, double massShift, boolean isStatic) {
        this.symbol = symbol;
        this.massShift = massShift;
        this.isStatic = isStatic;
    }

    public TerminalModification(char symbol, double massShift) {
        this.symbol = symbol;
        this.massShift = massShift;
    }
    public double getMassShift() {
        return massShift;
    }
    public char getSymbol() {
        return symbol;
    }
    public boolean isStatic() {
        return isStatic;
    }


}



