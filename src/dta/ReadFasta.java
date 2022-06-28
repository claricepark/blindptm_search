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
        FastaProtein protein = new FastaProtein();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        StringBuffer sb = new StringBuffer();

        String define = br.readLine(); //> line
        String description = define.substring(1);
        // if the line read is a empty string, ignore it
        while ((define = br.readLine()) != null && (define.equals("") || define.charAt(0) != '<')) {
            if(define.equals("")){
                continue;
            }else if (!define.equals("")) {
                String line = define.trim();
                sb.append(line);
            }else if(description.substring(0,1).equals(">")){
                proteinList.add(protein.setDescription(sb));
            }else{
                proteinList.add(protein.setSequence(sb));
            }
        }

        System.out.println(proteinList);
        sb.setLength(0);
        br.close();
    }

}

