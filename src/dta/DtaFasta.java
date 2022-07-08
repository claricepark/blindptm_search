package dta;

import java.io.*;
import java.util.*;


public class DtaFasta {

    public static void main(String[] args) throws Exception { //cannot be main method
    }
    public static String getIdenfitiedFasta(String dtaFile, String fastaFile) throws Exception {

//        Scanner inScanner = new Scanner(System.in);
//        System.out.print("Enter file path:"); ///Users/claricepark/data/blindptm/search.xml
//        String DtaPathInput = inScanner.next();
//        String inputDtaFile = DtaPathInput; //"/Users/claricepark/data/blindptm/DTASelect-filter.txt";
//        System.out.println("You entered: " + DtaPathInput);
        Set<String> dtaProteinList = DtaReader.getProteinList(dtaFile);

//        Scanner sc = new Scanner(System.in);
//        System.out.print("Enter second file path:");
//        String fastaPathInput = sc.next();
//        String inputFastaFile = fastaPathInput;//"/Users/claricepark/data/blindptm/UniProt_human_reviewed_contaminant_05-23-2020_reversed.fasta";
//        System.out.println("You entered: " + fastaPathInput);
        List<FastaProtein> fastaProteinList = ReadFasta.getFastaProteinList(fastaFile);

        String outputPath = fastaFile.substring(0, fastaFile.lastIndexOf('/')) + File.separator + "identified.fasta";
        //subset fasta file: "/Users/claricepark/data/blindptm/identified_subset.fasta";
        //String inputFasta = "/Users/claricepark/data/blindptm/small.fasta";

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
        System.out.println("Output file successfully created");

        for(String protein : dtaProteinList){

            if(fastaProteinSet.contains(protein)) {
                FastaProtein fp = fastaHt.get(protein);
                //System.out.println("Output file successfully created");
//                System.out.println(fp.getDescription());
//                System.out.println(fp.getSequence());
                writer.write(">" + fp.getDescription());
                writer.write(System.getProperty( "line.separator" ));
                writer.write(fp.getSequence());
                writer.write(System.getProperty( "line.separator" ));
            }
        }

        writer.close();
        return("Output file successfully created");
    }
}
