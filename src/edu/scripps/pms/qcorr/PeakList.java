/**
 * @file PeakList.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeakList
 * @author Tao Xu
 * @date $Date: 2007/02/27 01:31:02 $
 */



package  edu.scripps.pms.qcorr;

import java.util.ArrayList;
import java.util.ListIterator;
import java.util.Iterator;
import java.util.List;
import java.util.Comparator;
import java.util.Collections;

public class PeakList {

    public static int DEFAULTNUMPEAKS = 1000;
    public static final int DEFAULTSPECTRUMSIZE = 100000;   
    public static boolean SORTBYM2Z = false;
    public static boolean SORTBYINTENSITY = true;
    private float totalIntensity = 0;
    private int loscan;
    private int hiscan;
    private String scannum;
    private float precursorMass;
    private List<String> hlines;
    private List<String> ilines = new ArrayList<String>();
    private ArrayList<Zline> zlines = new ArrayList<Zline>();
    private ArrayList<Peak> peaks = new ArrayList<Peak>(DEFAULTNUMPEAKS);
    private String listType;
    private int index = 0;

    private ArrayList<Peak> peaksSortedByIntensity = null;
    private ArrayList<Peak> peaksSortedByM2z = null;

    public int getNumZlines() {
        return zlines.size();
    }
    public void setHlines(List<String> hlines) {
       this.hlines = hlines;
    }
    public void setZlines(ArrayList<Zline> zlines) {
       this.zlines = zlines;
    }
    public void setLoscan(int loscan) {
        this.loscan = loscan;
    }
    public void setHiscan(int hiscan) {
        this.hiscan = hiscan;
    }
    public List<String> getIlines() {
        return ilines;
    }
    public List<String> getHlines() {
        return hlines;
    }
    public void setPrecursorMass(float precursorMass) {
        this.precursorMass = precursorMass;
    }
    public void addIline(String l) {
        ilines.add(l);
    }
    public void addZline(Zline z) {
        zlines.add(z);
    }

    public void addPeak(Peak p) {
        p.setIndex(index++);
        peaks.add(p);
        totalIntensity += p.getIntensity();
    }

    public int numPeaks() {
        return peaks.size();
    }
    public float getTotalIntensity() {
        return totalIntensity;
    }
    public Iterator<Zline> getZlines() {
        return zlines.listIterator();
    }
    public void sortPeaks(boolean sortByIntensity) {
        peaks = getSortedPeaks(sortByIntensity);
    }

    /**
     * Return a sorted list (incremental by m2z or intensity)
     * for the Peaks in this PeakList
     * @param sortByIntensity - indicate how the list should be sorted,
     *                          true for sort by intensity,
     *                          false for sort by M2z
     * @note user must not modify the List returned. 
     * 
     */
    public synchronized ArrayList<Peak> getSortedPeaks(boolean sortByIntensity) {
        ArrayList<Peak> sortedPeaks =  new ArrayList(DEFAULTNUMPEAKS);
        if (sortByIntensity) {
            if (peaksSortedByIntensity != null) {
                return peaksSortedByIntensity;
            } else {
                peaksSortedByIntensity = sortedPeaks;
            }
        } else {
            if (peaksSortedByM2z != null) {
                return peaksSortedByM2z;
            } else {
                peaksSortedByM2z = sortedPeaks;
            }
        }
        for (Peak p : peaks) {
            sortedPeaks.add(p);
        }

        Collections.sort(sortedPeaks, new PeakComparator(sortByIntensity));
        return sortedPeaks; 
    }
    public synchronized ArrayList<Peak> getSortedPeaks(int numPeaks, boolean sortByIntensity) {
        ArrayList topPeaks = new ArrayList<Peak>(numPeaks);
 
        List<Peak> sortedList = getSortedPeaks(sortByIntensity);
        int totalPeaks = sortedList.size()-1; 
        while(numPeaks > 0 && totalPeaks >= 0) {
            Peak p = sortedList.get(totalPeaks--);
            //if(p.getM2z() > getMaxChargeState()*getPrecursorMass()){
            topPeaks.add(p);
            numPeaks--;
             
        }
        return topPeaks;
    }
    public ListIterator<Peak> getPeaks() {
        return peaks.listIterator();
    }

    public String getListType()
    {
        return listType;
    }

    public void setListType(String listType)
    {
        this.listType = listType;
    }
    /**
     * Return the low scan number of this
     */
    public int getLoscan() {
        return loscan;
    }

    public int getHiscan() {
        return hiscan;
    }

    public float getPrecursorMass() {
        return precursorMass;
    }
    
    public float getMaxM2z() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        return sortedPeaks.get(numPeaks()-1).getM2z();
    } 
    public float getMinM2z() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        return sortedPeaks.get(0).getM2z();
    } 
    public float getMaxIntensity() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
        return sortedPeaks.get(sortedPeaks.size()-1).getIntensity();
    }
    public float getMinIntensity() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
        return sortedPeaks.get(0).getIntensity();
    }
    public String getSpectrumWithHlines() {
        StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
        for(String h : hlines) {
            sb.append(h);
            sb.append("\n");
        }
        getSpectrumWithoutHlines(sb); 
        return sb.toString();       
    }
    public String getSpectrumWithoutHlines() {
        StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
        getSpectrumWithoutHlines(sb);
        return sb.toString();        
    }
    private void getSpectrumWithoutHlines(StringBuffer sb) {
        if(loscan<10){
	    scannum = "00000" + new Integer(loscan).toString();
	} else if (loscan>=10 && loscan<100 ) {
	    scannum = "0000" + new Integer(loscan).toString();	
	} else if (loscan>=100 && loscan<1000 ) {
	    scannum = "000" + new Integer(loscan).toString();
	} else if (loscan>=1000 && loscan<10000 ) {
	    scannum = "00" + new Integer(loscan).toString();
	} else if (loscan>=10000 && loscan<100000 ) {
	    scannum = "0" + new Integer(loscan).toString();
	} else {
	    scannum = new Integer(loscan).toString();
	}
	addSline(sb);
        addIlines(sb);
        addZlines(sb);
        addSpectrum(sb);
    }
    private void addIlines(StringBuffer sb) {
        for(String i : ilines) {
            sb.append(i);
            sb.append("\n");
        }
    }
    private void addSline(StringBuffer sb) {
        sb.append("S");
        sb.append("\t");
        sb.append(scannum);
        sb.append("\t");
        sb.append(scannum);
        sb.append("\t");
        sb.append(precursorMass);
        sb.append("\n");
    }
    private void addZlines(StringBuffer sb) {
        for(Zline z : zlines) {
            sb.append("Z");
            sb.append("\t");
            sb.append(z.getChargeState());
            sb.append("\t");
            sb.append(z.getM2z());
            sb.append("\n");
            for(String d : z.getDlines()) {
                sb.append(d);
                sb.append("\n");
            }
        }
    }
    private void addSpectrum(StringBuffer sb) {
        for(Peak p : peaks) {
            sb.append(p.getM2z());
            sb.append(" ");
            sb.append(p.getIntensity());
            sb.append("\n");
        }        
    }
    public int getMaxChargeState() {
        int maxChargeState = 0;
        for(Zline z : zlines) {
            if(z.getChargeState() > maxChargeState) {
                maxChargeState = z.getChargeState();
            }
        }
        return maxChargeState;
    } 
    private List<Peak> getMostIntensPeaks(int numPeaks, double lowM2zLimit, double highM2zLimit) {
        ArrayList<Peak> topPeaks = new ArrayList<Peak>(numPeaks);
 
        List<Peak> sortedList = getSortedPeaks(true);
        int totalPeaks = sortedList.size()-1; 
        while(numPeaks > 0 && totalPeaks >= 0) {
            Peak p = sortedList.get(totalPeaks--);
            //if(p.getM2z() > getMaxChargeState()*getPrecursorMass()){
            double m2z = p.getM2z();
            if(m2z > lowM2zLimit && m2z < highM2zLimit){
                topPeaks.add(p);
                numPeaks--;
            }
             
        }
        return topPeaks;
    }
    public PointList [] calcQCorrs(int numPeaks, int accuracyFactor, double massWindow) {
//        int massH = (int)(MassSpecConstants.MASSH*accuracyFactor + 0.5);
        PointList [] points = new PointList[getNumZlines()];
        int zlineCounter = 0;
        
        for(Zline z : zlines) {
            int chargeState = z.getChargeState();
            double lowLimit = massWindow*chargeState;
            double highLimit = precursorMass*chargeState - lowLimit;
            List<Peak> topList = getMostIntensPeaks(numPeaks, lowLimit, highLimit);
            int numBins = (int)(precursorMass*chargeState*accuracyFactor);
            double [] intens = new double[numBins];    
            for(Peak p : topList) {
                //int index = (int)(p.getM2z()*accuracyFactor+0.5);
                int index = (int)(p.getM2z()*accuracyFactor);
                intens[index] += p.getIntensity();
            } 
            //double [] temp = intens;
            intens = reverse(intens); // get the revered list for correlation
            int numShift = (int)(massWindow*chargeState*accuracyFactor);
            //int offSet = numShift/2; 
//System.out.println("scan#: " + getLoscan() +"\tprecuMass: " + precursorMass +  "\toffSet: " + offSet);
            int offSet = numShift/2+1;
	    numShift += 1;
            double [] qcorr = new double [numShift];
        
            for(Peak p : topList) {
                int m2z = (int)(p.getM2z()*accuracyFactor);
                double intensity = p.getIntensity();  
                int startIndex = m2z-offSet;
                for(int i = 0; i < numShift; i++) {
                    qcorr[i] += intensity*intens[startIndex + i];
                }
            } 
            PointList pointList = new PointList();
            double firstMass = precursorMass*chargeState+(massWindow/2*chargeState);
            double dAccFactor = accuracyFactor; // convert to accuracyFactor to double
//System.out.println("Max: " + (numBins+offSet)/dAccFactor + "\tMin: " + (numBins-offSet)/dAccFactor);
            for(int i = 0; i < qcorr.length; i++) {
                pointList.addPoint(new Point(firstMass-(i/dAccFactor), qcorr[i]));
            }
            points[zlineCounter++] = pointList;
            //int indexMax = getMaxIndex(qcorr);
            //double mass = (numBins-offSet-indexMax-1.5*massH)/(accuracyFactor+0.0); 
            //double mass = (firstMass-indexMax)/(accuracyFactor+0.0); 
//System.out.println("chargeState: " + z.getChargeState()+"\tnumBins: " + numBins +"\tm2z: " + z.getM2z() + "\tmass: " + mass + "\tindexMax: " + indexMax + "\tqcorr: " + qcorr[indexMax]);
//System.out.println(z.getChargeState() + "\t"+qcorr[indexMax] + "\t" + (z.getM2z()-mass));
            //String dline = "D\t" + mass + "\t" + qcorr[indexMax];
            //z.addDline(dline);
        } 
        return points;
    }
    private double [] reverse(double [] intens) {
        double [] reversed = new double[intens.length];
        int lastIndex = intens.length - 1;
        for(int i = 0; i <= lastIndex; i++) {
            reversed[lastIndex-i] = intens[i];
        }
        return reversed;
    }
    private int getMaxIndex(double [] qcorr) {
        
        int indexMax = 0;
        double max = 0;
        for(int i = 0; i < qcorr.length; i++) {
            if(qcorr[i] > max) {
                max = qcorr[i];
                indexMax = i;
            }            

        }
//System.out.println("max: " + max);
        return indexMax;
    } 
}
