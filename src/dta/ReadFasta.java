package dta;

import java.io.FileReader;
import java.util.*;
import java.io.BufferedReader;

public class ReadFasta {

    static List<FastaProtein> proteinList;

    public static void main(String[] args) throws Exception {

        String inputFile = "/Users/claricepark/data/blindptm/identified.fasta";
        //read fasta file.
        //ecah potein, create FastaProtein object and add to List (ArrayList)

        //List<FastaProtein> proteinList = new ArrayList();
        proteinList = new ArrayList<FastaProtein>();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
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

            proteinList.add(protein);
            sb.delete(0, sb.length());

            if(eachLine == null)
                break;
        }

        //System.out.println(proteinList);
        br.close();
    }

    public static List<FastaProtein> getList() {

        return proteinList;

    }

}

