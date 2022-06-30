package dta;

import java.util.*;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


public class DtaFasta {

    public static void main(String[] args) throws Exception {

        List<String> combine = new ArrayList();

        BufferedWriter bw = null;

        ReadFasta.getList();
        DtaReader.getHashSet();

        File textFile = new File("FastaMatch.txt");
        try {
            bw = new BufferedWriter(new FileWriter(textFile));
        } catch (IOException e) {
            e.printStackTrace();
        }

        String line;
        //int listsz = ReadFasta.proteinList.size();
        for (FastaProtein i : ReadFasta.proteinList) {
            for (FastaProtein j : DtaReader.data) {
                if (i.equals(j)) {
                    combine.add(String.valueOf(i));

                }
                try {
                    if (bw != null)
                        bw.close();
                } catch (IOException e) {
                    //
                }
            }
        }System.out.println(combine);
        bw.write(combine + "\n");

    }

}
