package scripts.imscore;

import edu.scripps.pms.util.io.SQTParser;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.sqt.SQTPeptide;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.Random;

public class SQTTargetDecoyDeltaIMobilityDistribution {
    public static void main(String args[]) throws IOException
    {
//        scoreDistribution("/home/ip2/data/timstof_hela_yates/ty_asms2020/190806_200ng_180m_01/timscore/align/190806_200ng_180m_01_Slot2-3_1_632_nopd.sqt");
        scoreDistribution("/home/ip2/data/ccs_score_dist/test.sqt");


//        evaluate(args[0]);
    }


    public static void scoreDistribution(String filename) throws IOException {

        DecimalFormat df=new DecimalFormat("0.0");

        SQTParser reader = new SQTParser(filename);

        SQTPeptide peptide;
        MLine mLine;

        int count=0;

        for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
        {
            peptide = itr.next();

            String[] arr = peptide.getSLine().split("\t");
            double expIMValue = Double.parseDouble(arr[arr.length-1]);

            StringBuffer sb = new StringBuffer();
            double xcorr = 0;
            double deltaCN = 10;
            String topProtein = "";
            String secondProtein = "";
            double timScore = 0;

            for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
            {
                mLine = mItr.next();
                String mLineString = mLine.getMLine();
                String[] mArr = mLineString.split("\t");

                if(!mLineString.startsWith("M\t1"))
                    continue;
                String timScoreStr = mArr[mArr.length-1];
                if(timScoreStr.equals("NA"))
                    continue;

//                if(mArr.length<13 || mArr[mArr.length-1].equals("NA"))
//                    continue;

                xcorr = mLine.getXcorrValue();
                if(xcorr<=0)
                    continue;

                deltaCN = mLine.getDeltaCnValue();
                timScore = Double.parseDouble(timScoreStr);


                StringBuffer proteinSb = new StringBuffer();
                boolean isTarget=false;


                for(Iterator<String> lItr = mLine.getLLine(); lItr.hasNext(); )
                {
                    String protein = lItr.next();
                    proteinSb.append(protein).append(";");
                    if(!protein.startsWith("Rever"))
                        isTarget=true;
                }



                if(!isTarget && timScore > 0.7) {
                    Random random = new Random();
                    double value = random.nextDouble();
                    if(value>0.1)
                        continue;
                }

                //System.out.println( (isTarget?1:0) + "\t" +  df.format(xcorr) + "\t" + df.format(deltaCN) + "\t" +  df.format(timScore) );
                System.out.println( (isTarget?1:0) + "\t" + xcorr + "\t" + deltaCN + "\t" +  timScore );
            }
        }
    }


    /**********************
     * target and decoy distribution
     */
    public static void evaluate(String filename) throws IOException {

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

                if(!mLineString.startsWith("M\t1"))
                    continue;


                if(mArr.length<13 || mArr[mArr.length-1].equals("NA"))
                    continue;


                double topMValue = Double.parseDouble(mArr[mArr.length-2]);

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

                double deltaIm = (topMValue - expIMValue)/((topIMValue>expIMValue)?topMValue:expIMValue);
                deltaIm = Math.abs(deltaIm);
                System.out.println(isTarget + "\t" +  deltaIm + "\t" + proteinSb.toString());


            }

                count++;
            }
        }



}
