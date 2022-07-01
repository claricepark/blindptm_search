package dta;

import java.util.*;
import java.io.BufferedReader;
import java.io.FileReader;

public class DtaReader {

    public static Set<String> getProteinList(String dtaFile) throws Exception {

        Set<String> proteinList = new HashSet<String>();

       // String inputFile = "/Users/claricepark/data/blindptm/DTASelect-filter.txt";
        //String inputFile = args[0];
        String eachLine = null;

        BufferedReader br = new BufferedReader(new FileReader(dtaFile));
        while ((eachLine = br.readLine()) != null) {
            String[] arr = eachLine.split("\t");


            if(arr.length >= 8 && arr.length <= 15){
                if(!arr[0].equals("Locus")) {
                    //System.out.println(eachLine);
                    proteinList.add(arr[0]);
                }
            }
        }
//        System.out.println(data);

        br.close();

        return proteinList;
    }

}

