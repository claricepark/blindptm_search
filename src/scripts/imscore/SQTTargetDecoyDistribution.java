package scripts.imscore;

import edu.scripps.pms.util.io.SQTParser;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.sqt.SQTPeptide;

import java.io.IOException;
import java.util.Iterator;

public class SQTTargetDecoyDistribution {
    public static void main(String args[]) throws IOException
    {
        evaluate(args[0]);
    }

    /**********************
     * target and decoy distribution
     */
    public static void evaluate(String filename) throws IOException {

        //SQTParser reader = new SQTParser( "/home/ip2/data/hela_qc/ms2/small.sqt" );
        // SQTParser reader = new SQTParser( "/home/ip2/data/hela_qc/ms2/hela_qc_nopd.sqt" );
        System.out.println("file name\t" + filename);

        SQTParser reader = new SQTParser(filename);

        SQTPeptide peptide;
        MLine mLine;

        System.out.println("start...");

        int count=0;

        for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
        {
            peptide = itr.next();

            String[] arr = peptide.getSLine().split("\t");
            double expIMValue = Double.parseDouble(arr[arr.length-1]);

            StringBuffer sb = new StringBuffer();

            double topXcorr = 0;
            double secondXcorr = 0;
            double deltaCN = 10;
            String topProtein = "";
            String secondProtein = "";
            double topIMScore = 0;
            double secondIMScore = 0;
            double topIMValue = 0;
            double secondIMValue = 0;

            for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
            {
                mLine = mItr.next();
                String mLineString = mLine.getMLine();
                String[] mArr = mLineString.split("\t");

                if(mArr.length<13 || mArr[mArr.length-1].equals("NA"))
                    continue;

                StringBuffer proteinSb = new StringBuffer();
                boolean isTarget=false;
                for(Iterator<String> lItr = mLine.getLLine(); lItr.hasNext(); )
                {
                    String protein = lItr.next();
                    proteinSb.append(protein).append(";");
                    if(!protein.startsWith("Rever"))
                        isTarget=true;
                   //     secondProtein += temp + ";";

                 //   sb.append(temp + ";"); //System.out.print(temp + ";");
                }
                //System.out.println(isTarget + "\t" +  mArr[mArr.length-2] + "\t" + mArr[mArr.length-1] + "\t" + proteinSb.toString());
                System.out.println(isTarget + "\t" + mArr[mArr.length-1] + "\t" + proteinSb.toString());
            }
                ///System.out.println(deltaCN +  "\timvalue (" +   expIMValue + "\t"+ topIMValue + "\t" + secondIMValue + ")\t(" + topXcorr + "\t" + secondXcorr + ")\t(" + topIMScore + "\t" + secondIMScore + ")\t" + sb.toString());
               // System.out.println(expIMValue + "\t"+ topIMValue + "\t" + secondIMValue + "\t" + topXcorr + "\t" + secondXcorr + "\t" + sb.toString());
                count++;
            }
        }

//        System.out.println("total count\t" + count);



}
