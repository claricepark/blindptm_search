
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2008</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version $Id
 */

package edu.scripps.pms.util.io;

import edu.scripps.pms.util.io.SQTParser;
import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.seq.Fasta;

public class SqtReader {
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;


    
    public SqtReader(String fileName) throws IOException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        br = new BufferedReader(new FileReader(fileName));

        readHeader();
    }


    private void readHeader() throws IOException
    {

        while( (lastLine = br.readLine()).startsWith("H\t") );  //read head line
        while( (lastLine = br.readLine()).equals("") ); //read empty line
    }

    public Iterator<SQTPeptide> getSQTPeptide() throws IOException
    {
        return new Iterator<SQTPeptide> ()
        {
            private SQTPeptide peptide;
            private MLine mLine;

            public boolean hasNext()
            {
                return lastLine != null;
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

                peptide = new SQTPeptide(lastLine);


                java.util.ArrayList<MLine> list = new java.util.ArrayList<MLine>();
                
                while ( ( (lastLine = br.readLine()) != null && !lastLine.equals("")) &&
                        (lastLine.startsWith("M\t") || lastLine.startsWith("L\t")) )
                {
                    //sqt file has wrong deltaCN value.  So, grep the right value from it.                
                    boolean isDeltaCN=false;
                    double deltaCN=0;                
                    double temp;

                    list.clear();
                    
                    while(lastLine!=null && lastLine.startsWith("M\t"))
                    {
			//System.out.println(lastLine);
			/*  If M line is messed up in sqt file, discard it */
			String[] arr = lastLine.split("\t");
			if(arr.length<10)
			{
				isValidPeptide = false;
				break;
			}

                        mLine = new MLine(arr);

                        if(!isDeltaCN)
                        {
                            temp = Double.parseDouble( mLine.getDeltCN() );                        
                        
                            if(temp>0)
                            {
                                deltaCN=temp;                            
                                isDeltaCN=true;
                            }
                        }                                                
                        
                        while ( (lastLine = br.readLine()) != null &&  lastLine.startsWith("L\t") )
                        {
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

                return peptide;
            }
        };
    }

    public void close() throws IOException {
        br.close();
    }

}
