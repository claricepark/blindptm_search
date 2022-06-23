/*
 * ModProtein.java
 *
 * Created on May 3, 2005, 3:39 PM
 */

package edu.scripps.pms.util.dtaselect;

import edu.scripps.pms.util.seq.Fasta;

/**
 *
 * @author rpark
 */
public class ModProtein {
    
    /** Creates a new instance of ModProtein */
    private String locus;
    
    public ModProtein(String proteinLine) throws ArrayIndexOutOfBoundsException
    {
        this( proteinLine.split("\t") );
    }

    public ModProtein(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        this.locus = Fasta.getAccession(strArr[1]);
        
        
        //peptideList = new ArrayList<Peptide>();

        //init(strArr);
    }

    public String getLocus() {
        return locus;
    }

    public void setLocus(String locus) {
        this.locus = locus;
    }
    
}
