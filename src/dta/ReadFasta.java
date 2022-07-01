package dta;

import java.io.FileReader;
import java.util.*;
import java.io.BufferedReader;

public class ReadFasta {

    public static void main(String[] args) throws Exception {
        getFastaProteinList( "/Users/claricepark/data/blindptm/small.fasta");
    }

    public static List<FastaProtein> getFastaProteinList(String fastaFile) throws Exception {

        //String inputFile = "/Users/claricepark/data/blindptm/small.fasta";

        //List<FastaProtein> proteinList = new ArrayList();
        List<FastaProtein> proteins = new ArrayList();
        BufferedReader br = new BufferedReader(new FileReader(fastaFile));
        StringBuffer sb = new StringBuffer();

        String eachLine = null;
        //1. create FastaProtein object
        //2. set values correctly
        //3. add them to proteinList

        eachLine = br.readLine();
        while (true) {

            String defLine = eachLine.substring(1);
            while(   ((eachLine = br.readLine()) != null) && !eachLine.startsWith(">") ) {
                if(eachLine.equals(""))
                    continue;
                sb.append(eachLine);
            }
            FastaProtein protein = new FastaProtein();
            protein.setDescription(defLine);
            protein.setSequence(sb.toString());

            proteins.add(protein);
            sb.delete(0, sb.length());

            if(eachLine == null)
                break;
        }

        //System.out.println(proteinList);
        br.close();
        return proteins;
    }



}

