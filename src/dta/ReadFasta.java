package dta;

import org.apache.axis2.context.externalize.DebugObjectInput;

import java.io.FileReader;
import java.util.*;
import java.io.BufferedReader;

public class ReadFasta {

    private String lastLine;

    public static void main(String[] args) throws Exception {

        String inputFile = "/Users/claricepark/data/blindptm/identified.fasta";
        //read fasta file.
        //ecah potein, create FastaProtein object and add to List (ArrayList)

        List<ReadFasta> proteinList = new ArrayList();
        FastaProtein protein = null;
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        StringBuffer sb = new StringBuffer();

        String eachLine = null;
        //1. create FastaProtein object
        //2. set values correctly
        //3. add them to proteinList
        while ( (eachLine = br.readLine()) != null) {

            /*
            if(eachLine.equals(""))
                continue;

            String defLine = eachLine.substring(1);

            //while ( ((eachLine = br.readLine()) != null) && !eachLine.startsWith(">") )
            while(   ((eachLine = br.readLine()) != null) && !eachLine.startsWith(">") ) {

                sb.append(eachLine);
            }

            //new FastaProtein();
            */


            System.out.println(eachLine);
        }
        /*
        while ((define = br.readLine()) != null && (define.equals("") || define.charAt(0) != '<')) {
            if(define.equals("")){
                continue;
            }else if (!define.equals("")) {
                String line = define.trim();
                sb.append(line);
            }else if(description.substring(0,1).equals(">")){
                proteinList.add(protein.setDescription(sb.toString()));
            }else{
                proteinList.add(protein.setSequence(sb.toString()));
            }
        }
        */

        System.out.println(proteinList);
        br.close();
    }

}

