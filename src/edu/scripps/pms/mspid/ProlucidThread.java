/*
 *
 * @author Tao Xu
 * @email taoxu@scripps.edu
 * $Revision
 * $Date
 *
 */

package edu.scripps.pms.mspid;

import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.mspid.SearchResult;
import java.util.*;

public class ProlucidThread implements Runnable {
    private ProlucidSearchEngine pse; 
    private StringBuffer sb = new StringBuffer(20000);
    private int numprocessed = 0;
    public static final int NUMBUFFEREDRESULT = 20;

    public ProlucidThread(ProlucidSearchEngine pe) {
        pse = pe;
    }

    public void run() {
//      action.threadCallback();

        try {
            PeakList pl = null;
            while((pl = pse.getSpectrum()) != null) {
//System.out.println("Now start processing scan " + pl.getLoscan());
//int tempnumprev = pl.numPeaks();
//pl.processHcdSpectrum();
//System.out.println("With HCD processing, number of peaks before: " + tempnumprev + " and after: " +  pl.numPeaks());
                ArrayList<SearchResult> results = pse.search(pl);
                numprocessed++;
                
                outputSearchResults(results); 
                if(numprocessed == NUMBUFFEREDRESULT) {
                    // output result
                    pse.outputSearchResults(sb);     
                    numprocessed = 0;
                    sb = new StringBuffer(20000);
                } 
                
//System.out.println("Now finish processing scan " + pl.getLoscan());
            }         

            // in case no more spectra and numprocessed did not reach NUMBUFFEREDRESULT
            if(numprocessed != 0) {
                // output result 
                pse.outputSearchResults(sb);
            }
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    private void outputSearchResults(List<SearchResult> results) {

        for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {

            SearchResult rest = resultIt.next();
            if(rest != null) {
                if(rest.getPrintOnlySLine())
                {
                    //continue;
                    //System.out.println(rest.outputResults());
                }
                sb.append(rest.outputResults());
            }
        }

    }


}

