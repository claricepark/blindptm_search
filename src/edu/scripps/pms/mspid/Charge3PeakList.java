/**
 * @file Charge3PeakList.java
 * This is the source file for edu.scripps.pms.mspid.Charge3PeakList
 * @author Tao Xu
 * @date $Date: 2007/03/09 20:54:16 $
 */
package edu.scripps.pms.mspid;

import java.util.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.mspid.PeptideHit;
import edu.scripps.pms.util.seq.*;


public class Charge3PeakList extends ProcessedPeakList {
    public Charge3PeakList(PeakList peaks, Zline z, SearchParams sp, MassCalculator mc) {
        super(peaks, z, sp, mc);
    }
}



