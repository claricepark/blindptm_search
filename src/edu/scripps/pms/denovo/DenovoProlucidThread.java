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

public class DenovoProlucidThread implements Runnable {
    private DenovoProlucid dnp; 

    public DenovoProlucidThread(DenovoProlucid dp) {
        dnp = dp;
    }

    public void run() {
//      action.threadCallback();

        try {
            PeakList pl = null;
            while((pl = dnp.getSpectrum()) != null) {
                dnp.search(pl);
            }         
        } catch(Exception e) {
            e.printStackTrace();
        }
    }
}

