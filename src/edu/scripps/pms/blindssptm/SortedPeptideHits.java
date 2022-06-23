/**
 * @file SortedLinkedList.java
 * This is the source file for edu.scripps.pms.blindptm.SortedLinkedList
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.blindssptm;




import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;

public class SortedPeptideHits {
    private int maxNumElement;
    private LinkedList<PeptideHit> hits = new LinkedList<PeptideHit>(); 
 
    public SortedPeptideHits(int limit) {
        maxNumElement = limit;
    }
    
    public List<PeptideHit> getPeptideHits() {
        return hits;
    }
    private PeptideHit addPeptideHit(int index, Fasta f, int start, int end) {

        if(index >= maxNumElement) {
            return null;
        }
        PeptideHit p = new PeptideHit(f, start, end);
        hits.add(index, p);
        if(hits.size() > maxNumElement) {
            hits.removeLast();
        }
        return p;
    }
 
    public PeptideHit addPeptideHit(double prob, Fasta f, int start, int end) {
        int index = findIndex(prob, 0, hits.size()-1); 
        return addPeptideHit(index, f, start, end);
       
    }
    private int findIndex(double p, int low, int high) {
        if(hits.size() == 0) {
            return 0;
        }
        int mid = 0;
        while (high > low) {
            mid = low + (high - low) / 2;
            double midP = hits.get(mid).getProbability();
            if (p < midP)
                high = mid;
            else if (p > midP)
                low = mid;
            else
                return mid;
           if((high-low) == 1) {
               if(p <= hits.get(low).getProbability()) {
                   return low;
               }
               if(p > hits.get(high).getProbability()) {
                   return high + 1;
               } else {
                   return high;
               }
           }
        }
        return mid;
    }
    private static int findIndex(double [] values, double p, int low, int high) {
        if(values.length == 0) {
            return 0;
        }
        int mid = 0;
        while (high > low) {
            mid = low + (high - low) / 2;
            double midP = values[mid];
            if (p < midP)
                high = mid;
            else if (p > midP)
                low = mid;
            else
                return mid;
            if((high-low) == 1) {
                if(p <= values[low]) {
                    return low;
                } 
                if( p > values[high]) {
                    return high + 1; 
                } else {
                    return high;                
                }
            }
        }
        return mid;
    }
    public static void main(String args[]) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            System.out.println("Input number of elements, or exit to finish");
            
            String line = br.readLine();
            if("exit".equals(line) || "quit".equals(line)) {
                System.exit(0);
            }
            int numElements = Integer.parseInt(line);
            System.out.println();
         
            double [] d = new double[numElements];
        while(true) {

            System.out.print("Input value for search, e.g., 0.39, or exit to finish: ");
            
            line = br.readLine();
            if("exit".equals(line) || "quit".equals(line) || "end".equals(line)) {
                break;
            }
            double value = Double.parseDouble(line);        
            for(int i = 0; i < d.length; i++) {
                d[i] = 0.1 * i; 
                System.out.print(d[i] + "\t");
  
            }
            System.out.println();
            System.out.println("Index for " + value + ": " + findIndex(d, value, 0, d.length-1));
        } 

    }
}
