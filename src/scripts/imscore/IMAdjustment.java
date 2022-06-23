package scripts.imscore;

import edu.scripps.pms.util.io.DTASelectFilterReader;

import java.io.*;

public class IMAdjustment {
    public static void main(String args[]) throws IOException
    {
        String filename = args[0];
       // String path = filename.substring(filename.lastIndexOf(File.separator)) + File.separator;
        BufferedReader br = new BufferedReader(new FileReader(filename));
        //double imCorrectionValue = 0.02692342481708287;
//        double imCorrectionValue = -0.0269;
        double imCorrectionValue = Double.parseDouble(args[1]);

        String lastLine;
        PrintWriter out = new PrintWriter(filename + ".fixed");

        int count=0;
        while ((lastLine = br.readLine()) != null) {
            //System.out.println(lastLine);
//
//            if(count++>100)
//                break;

            if(!lastLine.startsWith("I\tIon M")) {
                out.println(lastLine);
                continue;
            }

            String[] arr = lastLine.split("\t");
            double imValue = Double.parseDouble(arr[2]);
            double correctedValue = imValue + imCorrectionValue;
            out.println(arr[0] + "\t" + arr[1] + "\t" + correctedValue);
//            System.out.println(lastLine + " " + correctedValue);

        }


        out.close();
        br.close();

    }
}
