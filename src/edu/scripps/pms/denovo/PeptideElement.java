/*
 *
 * @author Tao Xu
 * @email taoxu@scripps.edu
 * $Revision
 * $Date
 *
 */

package edu.scripps.pms.denovo;

import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.mspid.*;
import java.util.ArrayList;

public class PeptideElement {
    protected static SearchParams sp;
    protected String symbol;
    protected double mass;
    protected boolean isSimpleElement = true;
    public static final double MAXELEMENTMASS = 600;
    public static final int ACCURACYFACTOR = 1000; // the third digit
    public static final int MAXINDEX = (int)MAXELEMENTMASS * ACCURACYFACTOR;

    protected static ArrayList<ArrayList<PeptideElement>> mass2NTermElements = new ArrayList<ArrayList<PeptideElement>>(MAXINDEX);     
    protected static ArrayList<ArrayList<PeptideElement>> mass2CTermElements = new ArrayList<ArrayList<PeptideElement>>(MAXINDEX);     

    public PeptideElement(String symbol, double mass) {
        this.symbol = symbol;
        this.mass = mass;
    }

    public void addPeptideElement(PeptideElement pe) {
         
    }
    
    public ArrayList<PeptideElement> getNTermElements(double mass) {

        return mass2NTermElements.get((int)mass * ACCURACYFACTOR);
    }
    public ArrayList<PeptideElement> getCTermElements(double mass) {

        return mass2CTermElements.get((int)mass * ACCURACYFACTOR);
    }

    public static void setSearchParams(SearchParams p) {
        sp = p;
    }

    public double getMass() {
        return mass;
    } 

    public String getSymbol() {
        return symbol;
    }
}

