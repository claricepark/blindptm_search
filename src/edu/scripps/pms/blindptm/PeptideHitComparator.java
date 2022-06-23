
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2005</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu
 * @version 1.0
 * @date $Date: 2005/03/03 23:03:18 $
 */

package edu.scripps.pms.blindptm;

import java.util.Comparator;

public class PeptideHitComparator implements Comparator<PeptideHit> {

    public int compare(PeptideHit  ph1, PeptideHit ph2) {

        float f = ph1.getNumPeaksMatched() - ph2.getNumPeaksMatched();  // difference
         
        
        if (f > 0) {
            return 1;
        } else if (f < 0) {
            return -1;
        } else {
            return 0;
        }
    }

    /**
     * Return true if o is the same object as this, i.e., same memory address 
     */
    public boolean equals(Object o) {
        return this == o;
    }
}
