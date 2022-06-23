/*
 * BaseIdentificationReader.java
 *
 * Created on July 30, 2007, 1:23 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.util.io;

import java.io.IOException;

/**
 *
 * @author rpark
 */
public abstract class BaseIdentificationReader implements IdentificationReader {
    
    private String fileName;
        
    /** Creates a new instance of BaseIdentificationReader */
    public BaseIdentificationReader() {
    }
    
   
    
    
        //DTASelect specific methods
    public boolean isVersion2()
    {
        return false;
    }
    
    public double getConfidence()
    {
        return -1;
    }

/*
    public int getTotalPeptideNumber() throws IOException
    {
        int totalPeptideCount=0;

	Configuration conf = Configuration.getInstance();

	String idFileName = conf.getIdFileName();

	if(null == idFileName)
	    throw new Exception("Identification file is not found");

	if(conf.getIdFileName().endsWith(".txt"))
	{

	}
	else if(conf.getIdFileName().endsWith(".xml"))
	{
	    System.out.println("ee");
	}
        
        DTASelectFilterReader dtaReader = new DTASelectFilterReader(idFileName);

        for (Iterator<Protein> itr1 = dtaReader.getProteins(); itr1.hasNext(); ) {
            Protein protein = itr1.next();
            
            totalPeptideCount += protein.getPeptideSize();
        }
        
        return totalPeptideCount;


    }
    */
}
