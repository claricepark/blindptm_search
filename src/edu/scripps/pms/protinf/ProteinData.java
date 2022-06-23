/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.protinf;

import java.util.ArrayList;
import java.util.List;


/**
 *
 * @author harshil
 */
public class ProteinData 
{
    
    int groupId;
    int totalMembers;
    int occurence;
    private String sortingCode;
    int forwardProeinCount;
    int reverdeProteinCount;
    double globalRate;
    double confidence;
    int peptideItem;
    int totalSpecCount;
    double sumZScore;
    String uniquePeptide;
    String identifiedTrypticPeptide;
    String trypticPeptide;
    int proteinLength;
    double proteinSequenceCoverage;
    private boolean isReverseHit;
    String proteinAccession;
    String geneName;
    String defline;
    String allAccession;
    String allDefline;
    List<Integer> spectrumCount = new ArrayList();
    List<Integer> pepCount=new ArrayList();
    private List<Double> nsaf = new ArrayList();
    private List<Double> empai = new ArrayList();

    private List<Double> normSpcCount = new ArrayList();
    private List<Integer> searchIdGroup = new ArrayList<Integer>();
    private List<String> fileList = new ArrayList<String>();

    private double pValuetTest = 0;
    private double bHCorrection = 0;

    
   public List<String> getFileList()
	{
		return fileList;	
	}
    public void setFileList(List<String> file)
    {
	this.fileList = file;
    }
    public int getGroupId() {
        return groupId;
    }

    public void setGroupId(int groupId) {
        this.groupId = groupId;
    }

    public int getTotalMembers() {
        return totalMembers;
    }

    public void setTotalMembers(int totalMembers) {
        this.totalMembers = totalMembers;
    }

    public int getOccurence() {
        return occurence;
    }

    public void setOccurence(int occurence) {
        this.occurence = occurence;
    }

   

    public int getForwardProeinCount() {
        return forwardProeinCount;
    }

    public void setForwardProeinCount(int forwardProeinCount) {
        this.forwardProeinCount = forwardProeinCount;
    }

    public int getReverdeProteinCount() {
        return reverdeProteinCount;
    }

    public void setReverdeProteinCount(int reverdeProteinCount) {
        this.reverdeProteinCount = reverdeProteinCount;
    }

    public double getGlobalRate() {
        return globalRate;
    }

    public void setGlobalRate(double globalRate) {
        this.globalRate = globalRate;
    }

    public double getConfidence() {
        return confidence;
    }

    public void setConfidence(double confidence) {
        this.confidence = confidence;
    }

    public int getPeptideItem() {
        return peptideItem;
    }

    public void setPeptideItem(int peptideItem) {
        this.peptideItem = peptideItem;
    }

    public int getTotalSpecCount() {
        return totalSpecCount;
    }

    public void setTotalSpecCount(int totalSpecCount) {
        this.totalSpecCount = totalSpecCount;
    }

    public double getSumZScore() {
        return sumZScore;
    }

    public void setSumZScore(double sumZScore) {
        this.sumZScore = sumZScore;
    }

    public String getUniquePeptide() {
        return uniquePeptide;
    }

    public void setUniquePeptide(String uniquePeptide) {
        this.uniquePeptide = uniquePeptide;
    }

    public String getIdentifiedTrypticPeptide() {
        return identifiedTrypticPeptide;
    }

    public void setIdentifiedTrypticPeptide(String identifiedTrypticPeptide) {
        this.identifiedTrypticPeptide = identifiedTrypticPeptide;
    }

    public String getTrypticPeptide() {
        return trypticPeptide;
    }

    public void setTrypticPeptide(String trypticPeptide) {
        this.trypticPeptide = trypticPeptide;
    }

    public int getProteinLength() {
        return proteinLength;
    }

    public void setProteinLength(int proteinLength) {
        this.proteinLength = proteinLength;
    }

    public double getProteinSequenceCoverage() {
        return proteinSequenceCoverage;
    }

    public void setProteinSequenceCoverage(double proteinSequenceCoverage) {
        this.proteinSequenceCoverage = proteinSequenceCoverage;
    }

    public boolean isIsReverseHit() {
        return isReverseHit;
    }

    public void setIsReverseHit(boolean isReverseHit) {
        this.isReverseHit = isReverseHit;
    }

    public String getProteinAccession() {
        return proteinAccession;
    }

    public void setProteinAccession(String proteinAccession) {
        this.proteinAccession = proteinAccession;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getDefline() {
        return defline;
    }

    public void setDefline(String defline) {
        this.defline = defline;
    }

    public String getAllAccession() {
        return allAccession;
    }

    public void setAllAccession(String allAccession) {
        this.allAccession = allAccession;
    }

    public String getAllDefline() {
        return allDefline;
    }

    public void setAllDefline(String allDefline) {
        this.allDefline = allDefline;
    }

    public List<Integer> getSpectrumCount() {
        return spectrumCount;
    }

    public void setSpectrumCount(List spectrumCount) {
        this.spectrumCount = spectrumCount;
    }

    public void addSpectrumCount(String spectrumCount) {
        this.spectrumCount.add(Integer.parseInt(spectrumCount));
    }
    
    public List getPepCount() {
        return pepCount;
    }

    public void setPepCount(List pepCount) {
        this.pepCount = pepCount;
    }
    
    public void addPepCount(String pepCount) {
        this.pepCount.add(Integer.parseInt(pepCount));
    }

    public String getSortingCode() {
        return sortingCode;
    }

    public void setSortingCode(String sortingCode) {
        this.sortingCode = sortingCode;
    }

    public List<Double> getNsaf() {
        return nsaf;
    }

    public void setNsaf(List<Double> nsaf) {
        this.nsaf = nsaf;
    }

    public List<Double> getNormSpcCount() {
        return normSpcCount;
    }

    public void setNormSpcCount(List<Double> normSpcCount) {
        this.normSpcCount = normSpcCount;
    }

    public double getpValuetTest() {
        return pValuetTest;
    }

    public void setpValuetTest(double pValuetTest) {
        this.pValuetTest = pValuetTest;
    }

    /**
     * @return the bHCorrection
     */
    public double getbHCorrection() {
        return bHCorrection;
    }

    /**
     * @param bHCorrection the bHCorrection to set
     */
    public void setbHCorrection(double bHCorrection) {
        this.bHCorrection = bHCorrection;
    }

    /**
     * @return the searchIdGroup
     */
    public List<Integer> getSearchIdGroup() {
        return searchIdGroup;
    }

    /**
     * @param searchIdGroup the searchIdGroup to set
     */
    public void setSearchIdGroup(List<Integer> searchIdGroup) {
        this.searchIdGroup = searchIdGroup;
    }

    /**
     * @return the empai
     */
    public List<Double> getEmpai() {
        return empai;
    }

    /**
     * @param empai the empai to set
     */
    public void setEmpai(List<Double> empai) {
        this.empai = empai;
    }
    
    
    
    
}
