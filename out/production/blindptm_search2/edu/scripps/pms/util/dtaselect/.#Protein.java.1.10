package edu.scripps.pms.util.dtaselect;

/**
 * @author  Robin Park
 * @version $Id: Protein.java,v 1.10 2012/03/11 06:52:32 taoxu Exp $
 */
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

public class Protein
{
    private String locus;
    private String seqCount;
    private String spectrumCount;
    private String seqCoverage;
    private String length;
    private String molWt;
    private String pI;
    private String validation;
    private String description;
    private List<Peptide> peptideList;
    private String proteinLine="";
    private boolean isProblematic = false;

    public Protein(String proteinLine) throws ArrayIndexOutOfBoundsException
    {
        this( proteinLine.split("\t") );
        this.proteinLine = proteinLine;
    }

    public Protein(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        peptideList = new ArrayList<Peptide>();

        init(strArr);
    }

    public String getAccession() {
        return getLocus();
    }
    public Protein()
    {
        peptideList = new ArrayList<Peptide>();
    }
    public List<Peptide> getPeptideList() {
        return peptideList;
    }
    public void setElement(String[] strArr)
    {
	init(strArr);
        //System.out.println("in Protein, seqCoverage: " + seqCoverage);
    }

    public boolean addPeptide(Peptide peptide)
    {
        return peptideList.add(peptide);
    }

    public void removePeptide(Peptide peptide)
    {
        peptideList.remove(peptide);
    }

    public Iterator<Peptide> getPeptides()
    {
        return peptideList.iterator();
    }

    
    private void init(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        try {
            this.setLocus( Fasta.getAccession(strArr[0]) );
            this.setSeqCount(strArr[1]);
            this.setSpectrumCount(strArr[2]);
            this.setSeqCoverage(strArr[3]);
            this.setLength(strArr[4]);
            this.setMolWt(strArr[5]);
            this.setPI(strArr[6]);
            this.setValidation(strArr[7]);
            this.setDescription(strArr[8]);
        }
        catch(ArrayIndexOutOfBoundsException ex) {
            isProblematic = true;
            //throw new ArrayIndexOutOfBoundsException("Error : Mal formed DTASelect-filter.txt file in protein line: " + strArr[0]);
        }
    }
    public boolean isProblematic() {
        return isProblematic;
    }
    public String getLocus()
    {
        return locus;
    }

    public String getSeqCount()
    {
        return seqCount;
    }

    public String getSpectrumCount()
    {
        return spectrumCount;
    }

    public String getSeqCoverage()
    {
        return seqCoverage;
    }

    public String getLength()
    {
        return length;
    }

    public String getMolWt()
    {
        return molWt;
    }

    public String getPI()
    {
        return pI;
    }

    public String getValidation()
    {
        return validation;
    }

    public String getDescription()
    {
        return description;
    }

    public void setLocus(String locus)
    {
        this.locus = locus;
    }

    public void setSeqCount(String seqCount)
    {
        this.seqCount = seqCount;
    }

    public void setSpectrumCount(String spectrumCount)
    {
        this.spectrumCount = spectrumCount;
    }

    public void setSeqCoverage(String seqCoverage)
    {
        this.seqCoverage = seqCoverage;
    }

    public void setLength(String length)
    {
        this.length = length;
    }

    public void setMolWt(String molWt)
    {
        this.molWt = molWt;
    }

    public void setPI(String pI)
    {
        this.pI = pI;
    }

    public void setValidation(String validation)
    {
        this.validation = validation;
    }

    public void setDescription(String description)
    {
        this.description = description;
    }
    
    public int getPeptideSize()
    {
        return peptideList.size();
    }
    public int getNumPeptides()
    {
        return peptideList.size();
    }

    public String getProteinLine()
    {
        return proteinLine;
    }

}
