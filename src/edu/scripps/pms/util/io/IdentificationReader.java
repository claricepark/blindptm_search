/*
 * IdentificationReader.java
 *
 * Created on July 30, 2007, 11:45 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.util.io;

import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.dtaselect.Peptide;
import java.util.*;
import java.io.*;

/**
 *
 * @author rpark
 */
public interface IdentificationReader {

    //DTASelect specific methods
    public boolean isVersion2();
    public double getConfidence();
    public Iterator <Protein> getProteins() throws IOException;
    public int getTotalPeptideNumber() throws IOException;
    
    public String getFileName();
}
