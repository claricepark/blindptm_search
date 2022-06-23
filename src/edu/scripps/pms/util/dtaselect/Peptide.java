package edu.scripps.pms.util.dtaselect;

import java.util.ArrayList;

/**
 * @author  Robin Park, Tao Xu
 * @version $Id: Peptide.java,v 1.23 2014/08/12 19:13:59 rpark Exp $
 */
public class Peptide
{

    private boolean unique;
    private String fileName;
    private String xCorr;
    private String deltCN;
    private String mhPlus;
    private String calcMHplus;
    private String totalIntensity;
    private String spRank;
    private String spScore;
    private String ionProportion;
    private String redundancy;
    private String sequence;
    private String seqWithNoModification = null;
    private String[] peptideLine;
    private String tmpStr;
    private String scanNum;
    private String conf;
    private static int uniqueIndex = 0;
    private static int scanNumIndex = 1;
    private static int xcorrIndex = 2;
    private static int dcnIndex = 3;
    private static int confIndex = 4;
    private static int mPlusHIndex = 5;
    private static int calcMassIndex = 6;
    private static int totalIntensityIndex = 7;
    private static int spRankIndex = 8;
    private static int spScoreIndex = 9;
    private static int ionProportionIndex = 10;
    private static int redundancyIndex = 11;
    private static int sequenceIndex = 12;
    private static int pIIndex = 9; 
    private static int ppmIndex = -1; 
    private int hashcode = -1;
    private ArrayList<DiffModSite> diffmods = null;    
    private ArrayList<String> proteinAccessions = new ArrayList();
    private String filePath;

    public int hashCode() {
        if(hashcode == -1) {
            // fileName contains the scan number
            hashcode = (getSequence() + fileName).hashCode();
            //hashcode = (fileName + scanNum).hashCode();
            //System.out.println("file: " + fileName + "\tscannumber: " + scanNum);
        } 
//System.out.println(hashcode);
        return hashcode;
    } 
    public void setProteinAccessions(ArrayList<String> pacs) {
        proteinAccessions = pacs;
    }
    public ArrayList<String> getProteinAccessions() {
        return proteinAccessions;
    }
    public boolean equals(Object o) {
        Peptide p = (Peptide)o;
//System.out.println("in equals");
        return p == null? false : getSequence().equals(p.getSequence()) && fileName.equals(p.fileName);
    }

    public boolean isModifiedPeptide() {
        return sequence.indexOf("(") != -1 || sequence.indexOf("*") != -1 
                                           || sequence.indexOf("@") != -1
                                           || sequence.indexOf("#") != -1;
    }

    // get midSeq without diff mod info
    public String getSeqWithNoModification() {
         if(seqWithNoModification == null) {

            StringBuffer tempseq= new StringBuffer(getMidSeq());
            char c;
            double mod;

            //modifications
            double[] modifications= new double[tempseq.length()];
            for(int k=0; k<tempseq.length(); k++)
            {
                c= tempseq.charAt(k);
                if(c=='(')
                {
                    tempseq.delete(k, tempseq.indexOf(")")+1);
                    k--;
                }
                else if(!Character.isLetter(c))
                {
                        //mod=configuration[c];
                    tempseq.deleteCharAt(k);
                    k--;
                }
            }

            seqWithNoModification = tempseq.toString(); 

         }
         return seqWithNoModification;
    }
    public ArrayList<DiffModSite> getDiffModSites() {
        if (diffmods == null && isModifiedPeptide()) {

            StringBuffer tempseq= new StringBuffer(getMidSeq());
            diffmods = new ArrayList<DiffModSite>();
            char c;
            double mod;

            //modifications
            double[] modifications= new double[tempseq.length()];
            for(int k=0; k<tempseq.length(); k++)
            {
                c= tempseq.charAt(k);
                if(c=='(')
                {
                        mod=Double.parseDouble(tempseq.substring(k+1, tempseq.indexOf(")")));
                        if(k==0)
                            modifications[k]+=mod;
                        else
                            modifications[k-1]+=mod;
                        tempseq.delete(k, tempseq.indexOf(")")+1);
                        k--;
                }
                else if(!Character.isLetter(c))
                {
                        //mod=configuration[c];
                        mod=c;
                        if(k==0)
                            modifications[k]+=mod;
                        else
                            modifications[k-1]+=mod;
                        tempseq.deleteCharAt(k);
                        k--;
                }
            }

            seqWithNoModification = tempseq.toString(); 
            for(int i = 0; i < modifications.length; i++) {
                if(modifications[i] != 0){ 
                    //diffmods.add(new DiffModSite(i+1, modifications[i]));
                    diffmods.add(new DiffModSite(i+1, modifications[i], tempseq.charAt(i)));
                }
            }
        }
        return diffmods;
    }
    // this function need to be changed if the format of getInfo() function changes 
    public static String getInfoHeader() {
        return "Peptide\tXCorr\tDeltaCN\tMPlusH\tCalcMPlusH\tDeltaMass\tSpScore\tConfidence\tFileName";
    }
    public String getInfo() {
        StringBuffer sb = new StringBuffer(1000);
        sb.append(sequence + "\t" + xCorr + "\t" + deltCN);
        sb.append("\t" + mhPlus + "\t" + calcMHplus + "\t" + getDeltaMass());
        sb.append("\t" + spScore + "\t" + conf + "\t" + fileName);
        return sb.toString();
    }
    /* DTASelect 2.0 file */
    private void parseLine2() throws ArrayIndexOutOfBoundsException
    {
//System.out.println("confIndex: " + confIndex);
        scanNum = peptideLine[scanNumIndex];
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);
//System.out.println("scanNum: " + scanNum);
        this.setUnique((peptideLine[uniqueIndex]).startsWith("*"));
        this.setFileName(peptideLine[scanNumIndex]);
        this.setXCorr(peptideLine[xcorrIndex]);
        this.setDeltCN(peptideLine[dcnIndex]);
        this.setConf(peptideLine[confIndex]);
        this.setMhPlus(peptideLine[mPlusHIndex]);
        this.setCalcMHplus(peptideLine[calcMassIndex]);
        this.setTotalIntensity(peptideLine[totalIntensityIndex]);
        this.setSpRank(peptideLine[spRankIndex]);
        this.setSpScore(peptideLine[spScoreIndex]);
        this.setIonProportion(peptideLine[ionProportionIndex]);
        this.setRedundancy(peptideLine[redundancyIndex]);
        this.setSequence(peptideLine[sequenceIndex]);
//System.out.println("\t11: " + peptideLine[11] + "\t12: " + peptideLine[12]+ "\t13: " + peptideLine[13] + "\t14: " + peptideLine[14]);
    }
    
    public static void setFeatureIndices(String features) {
        String [] contents = features.split("\t");
        for(int i = 0; i < contents.length; i++) {
            String s = contents[i].trim();
            uniqueIndex = s.startsWith("Uni")? i : uniqueIndex; 
            scanNumIndex = s.startsWith("File")? i : scanNumIndex;
            xcorrIndex = s.startsWith("XC")? i : xcorrIndex;
            dcnIndex = s.startsWith("DeltCN")? i : dcnIndex;
            confIndex = s.startsWith("Conf%")? i : confIndex;
            mPlusHIndex = s.startsWith("M")? i : mPlusHIndex;
            calcMassIndex = s.startsWith("CalcM")? i : calcMassIndex;
            totalIntensityIndex = s.startsWith("Total")? i : totalIntensityIndex;
            spRankIndex = s.startsWith("SpR")? i : spRankIndex;
            spScoreIndex = s.startsWith("Prob")? i : spScoreIndex;
            ionProportionIndex = s.startsWith("IonP")? i : ionProportionIndex;
            redundancyIndex = s.startsWith("Red")? i : redundancyIndex;
            sequenceIndex = s.startsWith("Sequence")? i : sequenceIndex;
            pIIndex = s.startsWith("pI")? i : pIIndex;            
            ppmIndex = s.startsWith("PPM")? i : ppmIndex;            
            
        }
        //System.out.println("conf: " + confIndex);
    } 
    //For DTASelect version 2
    public Peptide(String[] peptideLine, boolean isV2)
    {
        this.peptideLine = peptideLine;
        
        //if(isV2)
            parseLine2();
        //else
            //parseLine();
    }

    
    /* DTASelect file */
    private void parseLine() throws ArrayIndexOutOfBoundsException
    {
        scanNum = peptideLine[1];
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);

        //this.setUnique(!"".equals(peptideLine[0]));
        this.setUnique((peptideLine[0]).startsWith("*"));
        this.setFileName(peptideLine[1]);
        this.setXCorr(peptideLine[2]);
        this.setDeltCN(peptideLine[3]);
        this.setMhPlus(peptideLine[4]);
        this.setCalcMHplus(peptideLine[5]);
        this.setTotalIntensity(peptideLine[6]);
        this.setSpRank(peptideLine[7]);
        this.setSpScore(peptideLine[8]);
        this.setIonProportion(peptideLine[9]);
        this.setRedundancy(peptideLine[10]);
        this.setSequence(peptideLine[11]);
    }
    public float getPi() {
        return Float.parseFloat(peptideLine[pIIndex]);
    }
    public float getDeltaMass() {
        
        return ppmIndex != -1 ? Float.parseFloat(peptideLine[ppmIndex]) : 1000;
    }
    public String getChargeState()
    {
        return this.fileName.substring( this.fileName.lastIndexOf(".") + 1 );
    }

    public String getLoScan()
    {
        tmpStr = this.fileName.substring( this.fileName.indexOf(".") +1 );

        return tmpStr.substring(0, tmpStr.indexOf(".") );
    }

    public String getFileName()
    {
        return fileName.substring(0, fileName.indexOf("."));
    }
    // return the fileName field that contains all the information
    public String getFileNameScanCharge() {
        return fileName;
    }

    public String getXCorr()
    {
        return xCorr;
    }

    public double  getXCorrValue()
    {
        return Double.parseDouble(xCorr);
    }
    public String getDeltCN()
    {
        return deltCN;
    }
    public double  getDeltCnValue()
    {
        return Double.parseDouble(deltCN);
    }

    public String getMhPlus()
    {
        return mhPlus;
    }

    public String getSpRank()
    {
        return spRank;
    }

    public String getSpScore()
    {
        return spScore;
    }
    public double getSpScoreValue()
    {
        return Double.parseDouble(spScore);
    }

    public String getIonProportion()
    {
        return ionProportion;
    }

    public String getRedundancy()
    {
        return redundancy;
    }

    public int getRedundancyValue()
    {
        return Integer.parseInt(redundancy);
    }
    public void increaseSpectrumCount()
    {
         int value = Integer.parseInt(redundancy);
         redundancy= String.valueOf(value+1);
                 
    }
    
    public String getSequence()
    {
        return sequence;
    }
    // return the peptide sequence without leading and tailing residues
    public String getMidSeq()
    {
        int lastindex = sequence.length() - 2;
        return sequence.substring(2, lastindex);
    }

    public boolean isUnique()
    {
        return unique;
    }

    public String getCalcMHplus()
    {
        return calcMHplus;
    }

    public String getTotalIntensity()
    {
        return totalIntensity;
    }

    public String getScanNum()
    {
        return scanNum;
    }
    public int getScanNumber()
    {
        return Integer.parseInt(scanNum);
    }

    public void setFileName(String fileName)
    {
        this.fileName = fileName;
    }

    public void setXCorr(String xCorr)
    {
        this.xCorr = xCorr;
    }

    public void setDeltCN(String deltCN)
    {
        this.deltCN = deltCN;
    }

    public void setMhPlus(String mhPlus)
    {
        this.mhPlus = mhPlus;
    }

    public void setSpRank(String spRank)
    {
        this.spRank = spRank;
    }

    public void setSpScore(String spScore)
    {
        this.spScore = spScore;
    }

    public void setIonProportion(String ionProportion)
    {
        this.ionProportion = ionProportion;
    }

    public void setRedundancy(String redundancy)
    {
        this.redundancy = redundancy;
    }

    public void setSequence(String sequence)
    {
        this.sequence = sequence;
    }

    public void setUnique(boolean unique)
    {
        this.unique = unique;
    }

    public void setCalcMHplus(String calcMHplus)
    {
        this.calcMHplus = calcMHplus;
    }

    public void setTotalIntensity(String totalIntensity)
    {
        this.totalIntensity = totalIntensity;
    }

    public void setScanNum(String scanNum)
    {
        this.scanNum = scanNum;
    }

    public String getConf() {
        return (null==conf)?"":conf;
    }
    public double  getConfValue() {
        if(conf == null) {
            return 0;
        }
        return Double.parseDouble(conf);
    }

    public void setConf(String conf) {
        this.conf = conf;
    }

    public String[] getPeptideLine()
    {
        return peptideLine;
    }
    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

}
