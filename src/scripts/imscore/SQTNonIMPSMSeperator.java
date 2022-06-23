package scripts.imscore;

import edu.scripps.pms.util.io.SQTParser;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.sqt.SQTPeptide;

import java.io.*;
import java.util.Iterator;

public class SQTNonIMPSMSeperator {
    public static void main(String args[]) throws IOException
    {
        runSeperator(args[0]);

        //System.out.println("aaa");
    }

    /**********************
     *
     */
    public static void runSeperator(String filename) throws IOException {

        StringBuffer headerSb = new StringBuffer();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String eachline;
        while( null != (eachline=br.readLine()) ) {
            if(!eachline.startsWith("H\t"))
                break;

            headerSb.append(eachline).append("\n");
            //System.out.println(eachline);
        }

        br.close();

        PrintWriter outIm = new PrintWriter(filename + ".im");
        PrintWriter outNoim = new PrintWriter(filename + ".noim");
        outIm.println(headerSb.toString());
        outNoim.println(headerSb.toString());

        System.out.println("file name\t" + filename);
        SQTParser reader = new SQTParser(filename);

        SQTPeptide peptide;
        MLine mLine;
        String temp;
        System.out.println("start...");

        int count=0;
        for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
        {
            peptide = itr.next();

            StringBuffer sb = new StringBuffer();
            sb.append(peptide.getSLine()).append("\n");

           // if(count++>10)
           //     break;

            if(peptide.getNumMlines()<=0)
                continue;

            String mLineStr = peptide.getMLine().next().getMLine();

            boolean imValid=false;

            if(mLineStr.endsWith("NA")) imValid=false;
            else imValid=true;

            //System.out.println(mLineStr);

            for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
            {

                mLine = mItr.next();
                sb.append(mLine.getMLine()).append("\n");

                for(Iterator<String> lItr=mLine.getLLine(); lItr.hasNext(); ) {
                    String lline = lItr.next();

                    sb.append("L\t").append(lline).append("\t0\t").append(mLine.getSequence()).append("\n");
                    //System.out.println(lline);
                }
            }

            if(imValid)
                outIm.println(sb.toString());
            else
                outNoim.println(sb.toString());
            //System.out.println(sb.toString());


        }
        outIm.close();
        outNoim.close();
    }


}
