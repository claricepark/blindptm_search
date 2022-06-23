package scripts.imscore;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

public class ConvertIMScoredSQTToOldSQT {
    public static void main(String args[]) throws IOException
    {
        String filename = args[0];
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String lastLine;
        PrintWriter out = new PrintWriter(filename + ".fixed");

        int count=0;
        while ((lastLine = br.readLine()) != null) {
       //     System.out.println(lastLine);

            if(!lastLine.startsWith("S\t") && !lastLine.startsWith("M\t")) {
                out.println(lastLine);
                continue;
            }

            String[] arr = lastLine.split("\t");

            if(lastLine.startsWith("S\t")) {
                for(int i=0;i<arr.length-1;i++) {
                    out.print(arr[i] + "\t");
                }
                out.println("");

            } else if(lastLine.startsWith("M\t")) {
                for(int i=0;i<arr.length-2;i++) {
                    out.print(arr[i] + "\t");
                }
                out.println("");

            } else {
                System.out.println("bug... fix me");
                System.out.println(lastLine);
            }


//            System.out.println(lastLine + " " + correctedValue);



        }


        out.close();
        br.close();

    }
}