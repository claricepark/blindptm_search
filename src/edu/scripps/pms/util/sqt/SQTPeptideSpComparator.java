
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
 * @date $Date: 2006/11/21 19:09:21 $
 */

package edu.scripps.pms.util.sqt;
import java.util.Comparator;

public class SQTPeptideSpComparator implements Comparator<SQTPeptide> {
 
    public static final int SORTBYXCORR = 1;
    public static final int SORTBYDELTACN = 2;
    public static final int SORTBYSP = 3;
    private int sortByScoreType;
 
    public SQTPeptideSpComparator() {
        this(SORTBYXCORR); // default sort by XCorr 
    }
    public SQTPeptideSpComparator(int scoreType) {
        sortByScoreType = scoreType;

    }
    public int compare(SQTPeptide s1, SQTPeptide s2) {
        double d = 0; 
        switch(sortByScoreType) {
            case 1 : d = s1.getTopHit().getXcorrValue() - s2.getTopHit().getXcorrValue(); break;
            case 2 : d = s1.getTopHit().getDeltaCnValue() - s2.getTopHit().getDeltaCnValue(); break;
            case 3 : d = s1.getTopHit().getSpValue() - s2.getTopHit().getSpValue(); break;
        }
        //jcom:w 
        
        if (d > 0) {
            return -1;
        } else if (d < 0) {
            return 1;
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
