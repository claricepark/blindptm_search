/**
 * @file Charge2PeakList.java
 * This is the source file for edu.scripps.pms.mspid.Charge2PeakList
 * @author Tao Xu
 * @date $Date: 2006/01/05 00:02:02 $
 */
package edu.scripps.pms.mspid;

import java.util.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.mspid.PeptideHit;
import edu.scripps.pms.util.seq.*;


class Charge2PeakList extends ProcessedPeakList {
    public Charge2PeakList(PeakList peaks, Zline z, SearchParams sp, MassCalculator mc) {
        super(peaks, z, sp, mc);
    }
}



