/**
 * @file Searchable.java
 * This is the source file for edu.scripps.pms.util.spectrum.Searchable interface
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.blindssptm;




import edu.scripps.pms.util.seq.Fasta;

import java.util.Iterator;
import java.io.IOException;

public interface Searchable {

    public Fasta accession2Fasta(String ac);
    public int getNumSequences(); 
    public Iterator<Fasta> getFastas();
    public Fasta getFasta(int index);
    
    public SearchResult search(ProcessedPeakList ppl) throws IOException;

}
