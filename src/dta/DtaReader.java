package dta;

import java.util.*;
import java.io.BufferedReader;
import java.io.FileReader;

public class DtaReader {
    public static void main(String[] args) throws Exception {

        String inputFile = "/Users/claricepark/data/blindptm/DTASelect-filter.txt";
        //String inputFile = args[0];
        String eachLine = null;
        Set<String> data = new HashSet<String>();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        while ((eachLine = br.readLine()) != null) {
            String[] arr = eachLine.split("\t");

            //remove locus
            //System.out.println(arr[0] + "\t" + arr.length);
            //System.out.println(eachLine);

            if(arr.length >= 8 && arr.length <= 13){
                if(!arr[0].equals("Locus")) {
                    //System.out.println(eachLine);
                    data.add(arr[0]);
                }
            }
        }
        System.out.println(data);

        br.close();
    }
//    public static Set<String> getHashSet() {
//
//        return data;
//    }
}

