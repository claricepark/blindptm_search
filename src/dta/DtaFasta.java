package dta;

import java.io.*;
import java.util.*;


public class DtaFasta {

    public static void main(String[] args) throws Exception {


        String inputFile = "/Users/claricepark/data/blindptm/DTASelect-filter.txt";
        String outputPath = inputFile.substring(0, inputFile.lastIndexOf('/')) + File.separator + "identified.fasta";

        //subset fasta file: "/Users/claricepark/data/blindptm/identified_subset.fasta";
        Set<String> dtaProteinList = DtaReader.getProteinList(inputFile);

        //String inputFasta = "/Users/claricepark/data/blindptm/small.fasta";
        String inputFasta = "/Users/claricepark/data/blindptm/UniProt_human_reviewed_contaminant_05-23-2020_reversed.fasta";
        List<FastaProtein> fastaProteinList = ReadFasta.getFastaProteinList(inputFasta);
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));


        Hashtable<String, FastaProtein> fastaHt = new Hashtable<>();
        for(FastaProtein fp:fastaProteinList) {
           // String desc = fp.getDescription().split(" ")[0];
            String key = "";
            if(fp.getDescription().length()>=21)
                key = fp.getDescription().substring(0,21);
            else
                key = fp.getDescription();

            fastaHt.put(key, fp);
        }

        Set<String> fastaProteinSet = fastaHt.keySet();

        FileWriter fw = new FileWriter(outputPath);

        for(String protein : dtaProteinList){

            if(fastaProteinSet.contains(protein)) {
                FastaProtein fp = fastaHt.get(protein);
                System.out.println(fp.getDescription());
                System.out.println(fp.getSequence());
                writer.write(fp.getDescription());
                writer.write(fp.getSequence());
            }
        }
    }

}
