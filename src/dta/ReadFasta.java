package dta;

import java.io.FileReader;
import java.util.*;
import java.io.BufferedReader;

public class ReadFasta {

    public static void main(String[] args) throws Exception {

        String inputFile = "/Users/claricepark/data/blindptm/identified.fasta";
        //read fasta file.
        //ecah potein, create FastaProtein object and add to List (ArrayList)

        List<ReadFasta> proteinList = new ArrayList();

        String eachLine = null;
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        while ((eachLine = br.readLine()) != null) {
            String[] arr = eachLine.split(">");

                FastaProtein protein = new FastaProtein();
//                protein.setDescription();
//                protein.setSequence();

                proteinList.add(protein.setDescription(arr[0]));


        }
        System.out.println(proteinList);
        br.close();
    }

}

