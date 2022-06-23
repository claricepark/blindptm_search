package edu.scripps.pms.protinf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.enzyme.Protease;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class ProteinTrypticPeptideNumber {



    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: proteintrypticpeptides <fastafile>");
	System.exit( -1);
    }



    public static int getNumTrypticPeptides(Fasta f, int minlength) {
        int counts = 0;
        int proteinlength =  f.getLength();
        int lastIndex = f.getLength() - 1;
        byte [] seq = f.getSequenceAsBytes();
            // mass of H2O need to be added
            // because the mass of the residues does not count H2O
        int templeft = 0;
        //int rightEnd = -1;
        int leftEnd = 0;

        while(leftEnd < proteinlength) {

            int rightEnd = leftEnd;
            
            while(rightEnd < lastIndex) {
                rightEnd++;
                if(rightEnd == lastIndex || seq[rightEnd-1] == 'K' ||  seq[rightEnd-1] == 'R') { // tryptic end
                    if(rightEnd - leftEnd > minlength-1) {
                        counts++;
                    }
                    break;
                }
            }
            leftEnd = rightEnd+1;

        }

        return counts;

    }
    //public static int getNumTrypticPeptides(Fasta f, Protease p, MassCalculator mc, double highLimit, double lowLimit, int minlength) {
    public static int getNumTrypticPeptides(Fasta f, MassCalculator mc, double highLimit, double lowLimit, int minlength) {
        double [] precMasses = mc.getPrecMasses();
        int counts = 0;
        int proteinlength =  f.getLength();
        int lastIndex = f.getLength() - 1;
        byte [] seq = f.getSequenceAsBytes();
            // mass of H2O need to be added
            // because the mass of the residues does not count H2O
        double mass = MassSpecConstants.MASSH3O;
        double tempmass = 0;
        int templeft = 0;
        //int rightEnd = -1;
        int leftEnd = 0;


//System.out.println("Protein Length: " + (lastIndex + 1) + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);




        while(leftEnd < proteinlength) {

            int rightEnd = leftEnd;
            
            while(rightEnd < lastIndex) {
                rightEnd++;
                if(rightEnd == lastIndex || seq[rightEnd-1] == 'K' ||  seq[rightEnd-1] == 'R') { // tryptic end
                    if(rightEnd - leftEnd > minlength-1) {
                        double precmass = mc.getPrecursorMass(seq, leftEnd, rightEnd);
                        if(precmass >= lowLimit && precmass <= highLimit) {
//System.out.println("Protein Length: " + (lastIndex + 1) + "\tleftEnd: " + leftEnd + "\trightEnd: " + rightEnd +"\tmass: " + precmass + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);

                            counts++;
                        }
                    }
                    break;
                }
            }
            leftEnd = rightEnd+1;

        }

        return counts;

    }

    public static void main(String[] args) throws Exception{

        
        HashSet<String> acset = new HashSet();
        ArrayList<String> aclist = new ArrayList();
        HashMap<String, Fasta> ac2protein = new HashMap();

        Protease p = new Protease("Trypsin");
        p.setType(true);

        p.addCleavageSite('R');
        p.addCleavageSite('K');

        MassCalculator mc = new MassCalculator();


	try {

            String fastaFileName = args[0];

           
             
            FileInputStream fis = new FileInputStream(new File(fastaFileName));
            for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
                Fasta f = it.next();
                String acc = f.getAccession();
                ac2protein.put(acc, f);
                aclist.add(acc);
                
            }
            fis.close();

	} catch(Exception e){
            help();
            e.printStackTrace(); 
	} finally {
	    //tidy up
	}
    }
}

