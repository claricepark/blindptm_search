package scripts.imscore;

import edu.scripps.pms.util.io.SQTParser;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.sqt.SQTPeptide;

import java.io.IOException;
import java.util.Iterator;

public class SQTEvaluator {
    public static void main(String args[]) throws IOException
    {
        filterPSMForGoodIMScore(args[0]);
    }

    /**********************
     * filter psm to detect good xcorr with low delta cn and good IMScore
     */
    public static void filterPSMForGoodIMScore(String filename) throws IOException {

        //SQTParser reader = new SQTParser( "/home/ip2/data/hela_qc/ms2/small.sqt" );
        // SQTParser reader = new SQTParser( "/home/ip2/data/hela_qc/ms2/hela_qc_nopd.sqt" );
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
                //System.out.println("seq===>>" + mLine.getSequence().replace('#', "");

                if(mLineString.startsWith("M\t1")) {
                    topXcorr = Double.parseDouble(mLine.getXcorr());
                    topIMScore  = Double.parseDouble(mArr[mArr.length-1]);
                    sb.append(peptide.getHiScan() + "\t" + peptide.getChargeState() + "\t" + expIMValue + "\t");
                    topIMValue =  Double.parseDouble(mArr[mArr.length-2]);
                }

                if(mLineString.startsWith("M\t2")) {
                    deltaCN = Double.parseDouble(mLine.getDeltCN());
                    secondXcorr = Double.parseDouble(mLine.getXcorr());
                    secondIMScore  = Double.parseDouble(mArr[mArr.length-1]);
                    secondIMValue =  Double.parseDouble(mArr[mArr.length-2]);
                }

//                    System.out.print(peptide.getHiScan() +  "\t" + peptide.getChargeState()  + "\t" + expIMValue + "\t");

                if(!mLineString.startsWith("M\t1") && !mLineString.startsWith("M\t2"))
                    continue;

                //System.out.println(mLine.getMLine());
                sb.append(mLine.getSequence() + "\t" + mLine.getXcorr() + "\t" + mArr[mArr.length-2] + "\t" + mArr[mArr.length-1] + "\t" );
                //System.out.print(mLine.getSequence() + "\t" + mLine.getXcorr() + "\t" + mArr[mArr.length-2] + "\t" + mArr[mArr.length-1] + "\t" );
                //System.out.println(mArr.length);

                for(Iterator<String> lItr = mLine.getLLine(); lItr.hasNext(); )
                {
                    temp = lItr.next();
                    if(mLineString.startsWith("M\t1"))
                        topProtein += temp + ";";
                    else if(mLineString.startsWith("M\t2"))
                        secondProtein += temp + ";";

                    sb.append(temp + ";"); //System.out.print(temp + ";");
                }

                sb.append("\t");
                //System.out.print("\t");

            }

            boolean isGoodScore=false;
            if(peptide.getChargeState().equals("2") && topXcorr>0.8)
                isGoodScore = true;
            else if(peptide.getChargeState().equals("3") && topXcorr>1.0)
                isGoodScore = true;

            if(isGoodScore && deltaCN<0.1) {
                //System.out.println(topProtein + " " + secondProtein);
                if(secondProtein.startsWith("Rever"))
                    continue;

                double im2ndDiff = Math.abs(expIMValue - secondIMValue);
                double imtopDiff = Math.abs(expIMValue - topIMValue);

                double imDiffAbs = Math.abs(secondIMValue - topIMValue);

                if(imDiffAbs<0.01)
                    continue;

                if(im2ndDiff> imtopDiff)
                    continue;

                //     if(secondIMScore>topIMScore  && (expIMValue - topIMValue)>0.02 && Math.abs(secondIMValue-topIMValue)>0.03)
                //if(secondIMScore<topIM && imDiff)
                //    continue;

                ///System.out.println(deltaCN +  "\timvalue (" +   expIMValue + "\t"+ topIMValue + "\t" + secondIMValue + ")\t(" + topXcorr + "\t" + secondXcorr + ")\t(" + topIMScore + "\t" + secondIMScore + ")\t" + sb.toString());
                System.out.println(expIMValue + "\t"+ topIMValue + "\t" + secondIMValue + "\t" + topXcorr + "\t" + secondXcorr + "\t" + sb.toString());
                count++;
            }
        }

        System.out.println("total count\t" + count);
    }


}
