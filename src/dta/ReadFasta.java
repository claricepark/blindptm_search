package dta;

import java.io.FileReader;
import java.util.*;
import java.io.BufferedReader;

public class ReadFasta {

    String inputFile = "/Users/claricepark/data/blindptm/identified.fasta";
    //read fasta file.
    //ecah potein, create FastaProtein object and add to List (ArrayList)

    List<ReadFasta> proteinList = new ArrayList();

    String eachLine = null;
    BufferedReader br = new BufferedReader(new FileReader(inputFile));
    while ((eachLine = br.readLine()) != null){
        ReadFasta seq = new ReadFasta();
        proteinList.add(seq);
    }
    //System.out.println(proteinList);


}