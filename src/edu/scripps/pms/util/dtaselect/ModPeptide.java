/*
 * ModPeptide.java
 *
 * Created on May 3, 2005, 4:02 PM
 */

package edu.scripps.pms.util.dtaselect;

/**
 *
 * @author rpark
 */
public class ModPeptide {
    
    private String modifiedResidue;
    private String massDifference;
    
    public ModPeptide(String peptideLine) throws ArrayIndexOutOfBoundsException
    {
        this( peptideLine.split("\t") );
    }

    public ModPeptide(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        //peptideList = new ArrayList<Peptide>();

        init(strArr);
    }
    
    private void init(String[] strArr)
    {
        this.modifiedResidue = strArr[1];
        this.massDifference = strArr[2];
        
        
    }

    public String getModifiedResidue() {
        return modifiedResidue;
    }

    public void setModifiedResidue(String modifiedResidue) {
        this.modifiedResidue = modifiedResidue;
    }

    public String getMassDifference() {
        return massDifference;
    }

    public void setMassDifference(String massDifference) {
        this.massDifference = massDifference;
    }
    
}
