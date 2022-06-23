package edu.scripps.pms.util.io;


import java.io.BufferedReader;
import java.util.Iterator;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.MLine;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version $Id: SQTParser.java,v 1.12 2006/01/18 00:53:46 taoxu Exp $
 */
public class SQTParser
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;

    public static void main(String args[]) throws IOException
    {
         SQTParser reader = new SQTParser( args[0] );

         SQTPeptide peptide;
         MLine mLine;
         String temp;
 	 System.out.println("start...");

         for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
         {
             peptide = itr.next();

             System.out.println("each====>>" + peptide.getSLine());


             for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
             {
                 mLine = mItr.next();

//                 System.out.println("seq===>>" + mLine.getSequence().replace('#', "");


                 for(Iterator<String> lItr = mLine.getLLine(); lItr.hasNext(); )
                 {
                     temp = lItr.next();
                     //System.out.println(temp);
                 }
             }
         }
    }

    public SQTParser(String fileName) throws IOException, NullPointerException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        init();
    }

    private void init() throws IOException, NullPointerException
    {
        br = new BufferedReader(new FileReader(fileName));

        readHeader();
    }

    private void readHeader() throws IOException, NullPointerException
    {
        lastLine = br.readLine();
        while(lastLine != null && !lastLine.startsWith("S")) {
            lastLine = br.readLine();
        }
/*
        if(lastLine.startsWith("H\t"))
        {
            while ((lastLine = br.readLine()).startsWith("H\t")); //read head line
            while ((lastLine = br.readLine()).equals("")); //read empty line
        }
        else
        {
            while (!(lastLine = br.readLine()).startsWith("S\t"));
        }
//System.out.println("lastLine: " + lastLine);
*/
    }

    public Iterator<SQTPeptide> getSQTPeptide() throws IOException
    {
        return new Iterator<SQTPeptide> ()
        {
            private SQTPeptide peptide;
            private MLine mLine;

            public boolean hasNext()
            {
                return lastLine != null && lastLine.startsWith("S");
            }

            public SQTPeptide next()
            {
                try
		{
                    peptide = getPeptide();
		    while(!isValidPeptide)
		    {
			isValidPeptide = true;
			peptide = getPeptide();
		    }
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return peptide;
            }

            public void remove()
            {
                throw new UnsupportedOperationException("Not supported");
            }

            private SQTPeptide getPeptide() throws IOException
            {

//                String[] strArr = lastLine.split("\t");

		//System.out.println("last line==>>" + lastLine);

		while( (lastLine != null) && lastLine.equals("") && !lastLine.startsWith("S\t") )
			lastLine = br.readLine();
//System.out.println(lastLine);
                peptide = new SQTPeptide(lastLine);


//System.out.println("Reached here. lastLine: " + lastLine);
                java.util.ArrayList<MLine> list = new java.util.ArrayList<MLine>();

                while ( ( (lastLine = br.readLine()) != null && !lastLine.startsWith("S")))
                {
                    //sqt file has wrong deltaCN value.  So, grep the right value from it.
//System.out.println("Reached here too. lastLine: " + lastLine);
                    boolean isDeltaCN=false;
                    float deltaCN=0;
                    float temp;

                    list.clear();

                    while(lastLine!=null && lastLine.startsWith("M\t"))
                    {
			//System.out.println(lastLine);
			/*  If M line is messed up in sqt file, discard it */
//System.out.println("lastLine: " + lastLine);
			String[] arr = lastLine.split("\t");
			if(arr.length<10)
			{
				isValidPeptide = false;
				break;
			}

                        mLine = new MLine(arr);
			mLine.setMLine(lastLine);

                        if(!isDeltaCN)
                        {
                            temp = Float.parseFloat( mLine.getDeltCN() );

                            if(temp>0)
                            {
                                deltaCN=temp;
                                isDeltaCN=true;
                            }
                        }

                        while ( (lastLine = br.readLine()) != null &&  lastLine.startsWith("L\t") )
                        {
			//System.out.println("---->>" + lastLine);
                            mLine.addLLine(lastLine);
                        }

                        list.add(mLine);
                    }

                    for(int i=0;i<list.size();i++)
                    {
                        list.get(i).setDeltCN(deltaCN);
                        peptide.addMLine(list.get(i));
                    }

                    //peptide.addMLine(mLine);
                    if(lastLine==null || lastLine.startsWith("S\t"))
                       break;
                }

////System.out.println("lastLine: " + lastLine);
                return peptide;
            }
        };
    }

    public void close() throws IOException
    {
        br.close();
    }

}
