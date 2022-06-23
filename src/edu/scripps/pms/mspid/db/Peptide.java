/**
 * @file Peptide.java
 * This is the source file for edu.scripps.pms.mspid.db.Peptide
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.mspid.db;



class Peptide implements Comparable<Peptide> {
    int mass;
    //byte [] location;
    int proteinIndex;
    char start;
    char stop;
    public Peptide(int mass, int proteinIndex, char start, char stop) {
        this.mass = mass;
        this.proteinIndex = proteinIndex;
        this.start = start;
        this.stop = stop;
    }
    public Peptide(int mass, byte[] loc) {
        this.mass = mass;
        //location = loc;
    }
    public int compareTo(Peptide p) {
        return mass - p.mass;
    }


}



