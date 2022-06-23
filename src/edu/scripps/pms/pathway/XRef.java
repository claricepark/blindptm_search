package edu.scripps.pms.pathway;

import java.io.BufferedReader;
import java.util.*;
import java.io.*;
import java.net.*;
//import edu.scripps.pms.util.sqt.SQTPeptide;
//import edu.scripps.pms.util.sqt.MLine;

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
 * @version $Id: XRef.java,v 1.3 2006/01/13 01:46:08 rpark Exp $
 */
public class XRef 
{
    private Hashtable<String, XProtein> proTable =null;
    private Hashtable<String, Set> pathTable = new Hashtable<String, Set>(); //Set : pathway name=protein set
    
    public List<XProtein> generatePathway(String xrefFile, List proteinList) throws IOException, Exception
    {
	proTable = MapFileGenerator.getXRef(xrefFile);

	XRef ref = new XRef();

        //List proteinList = new Vector();    
       
	//String lastLine;        

        for(Iterator<String> itr=proteinList.iterator(); itr.hasNext(); )
        {
            String each = itr.next();
	    //each = each.substring(0, each.indexOf('.'));
            //proteinList.add(each);
	    
	    XProtein protein = proTable.get(each);

	    if(null == protein)
	    {
		System.out.println(each + " is none");
		continue;
	    }

	    ref.readBioCarda( protein.getGeneList(), protein.getAccession() );
	}

          System.out.println(proteinList.size());
          System.out.println(proTable.size());
                
        List<XProtein> result = new Vector<XProtein>();
        
	for(Iterator<String> itr=proteinList.iterator(); itr.hasNext(); )
	{
            String accession = itr.next();            
	    XProtein protein = proTable.get(accession);

	    System.out.println("accession===>> " + accession + " " + protein);

	    //List l = protein.getGeneList();
            result.add(protein);
            
            /*
            System.out.println("accession==>>" + protein.getAccession());            
            
            for(Iterator<XGene> gitr=l.iterator(); gitr.hasNext(); )
	    {
		XGene eachGene = gitr.next();
		System.out.println("gene==>>" + eachGene.getName());            
                
                for(Iterator<Pathway> paitr=eachGene.getPathwayList().iterator(); paitr.hasNext(); )
                {
                    Pathway pway = paitr.next();

                    System.out.println("--->>" + pway.getName() + "\t" + pway.getUrl());
                }
            }
	      */    
	}
        
        return result;
        
    }
    
    public List generatePathway(String xrefFile, String inputFile) throws IOException, Exception
    {
        BufferedReader br = new BufferedReader(new FileReader(inputFile));

        String lastLine;        
        
        List<String> l = new Vector<String>();
        
	while( (lastLine=br.readLine())!= null )
	{
            String each = lastLine;

	    each = each.substring(0, each.indexOf('.'));
            l.add(each);            	    
	}        
        
        return generatePathway(xrefFile, l);
        //return l;        
    }
    
    public static void main(String args[]) throws IOException, Exception
    {
        XRef ref = new XRef();
        //ref.generatePathway(args[0], args[1]);
    }

    public static String readBioCarda(String geneName)
    {
	StringBuffer pathSb = new StringBuffer();

	try
	{
	    String url = "http://biocarta.com/genes/PathwayGeneSearch.asp?geneValue=" + geneName;

	    URL aURL = new URL(url);

	    URLConnection uc = aURL.openConnection();
	    BufferedReader in = new BufferedReader(
		    new InputStreamReader(
			uc.getInputStream()));
	    String inputLine;
	    StringBuffer sb = new StringBuffer();

	    while ((inputLine = in.readLine()) != null)
	    {
		sb.append(inputLine);
	    }

	    in.close();

	    if( sb.indexOf("No Pathway") >=0 )
		return "";	

	    String result = sb.toString();

	    while(true)
	    {
		int index = result.indexOf("bullet.gif");
		if(index<0)
		    break;

		result = result.substring( result.indexOf("bullet.gif")+6 );
		String pathUrl = result.substring( result.indexOf("href=\"")+6 );

		pathUrl = pathUrl.substring(0, pathUrl.indexOf("\">") );
		String pathway = result.substring( result.indexOf("href"), result.indexOf("</a>") );

		pathway = pathway.substring( pathway.indexOf('>')+1 );

		//            eachGene.addPathway(pathway, pathUrl);

		pathSb.append(pathway).append(",").append(pathUrl).append(";");
	    }
	}
	catch(Exception e)
	{
	    System.out.println("Error " + e);
	}

	return pathSb.toString();
    }

    public void readBioCarda(List geneList, String proteinAcc)
    {
	try
	{
	    for(Iterator<XGene> itr=geneList.iterator(); itr.hasNext(); )
	    {
		XGene eachGene = itr.next();
                
		String url = "http://biocarta.com/genes/PathwayGeneSearch.asp?geneValue=" + eachGene.getName();

		URL aURL = new URL(url);

		URLConnection uc = aURL.openConnection();
		BufferedReader in = new BufferedReader(
			new InputStreamReader(
			    uc.getInputStream()));
		String inputLine;
		StringBuffer sb = new StringBuffer();

		while ((inputLine = in.readLine()) != null)
		{
		    sb.append(inputLine);
		}

		in.close();

		if( sb.indexOf("No Pathway") >=0 )
		    continue;

		String result = sb.toString();

		while(true)
		{
		    int index = result.indexOf("bullet.gif");
		    if(index<0)
			break;

		    result = result.substring( result.indexOf("bullet.gif")+6 );
                    String pathUrl = result.substring( result.indexOf("href=\"")+6 );

                    pathUrl = pathUrl.substring(0, pathUrl.indexOf("\">") );
		    String pathway = result.substring( result.indexOf("href"), result.indexOf("</a>") );
                    
		    pathway = pathway.substring( pathway.indexOf('>')+1 );
		    System.out.println("2===>>" + pathway);
                    System.out.println("3===>>" + pathUrl);
                    
                    eachGene.addPathway(pathway, pathUrl);

		    Set tmpSet = pathTable.get(pathway);
		    if(null==tmpSet)
		    {
			tmpSet = new HashSet();
			tmpSet.add(proteinAcc);
			pathTable.put(pathway, tmpSet);
		    }
		    else
		    {
			tmpSet.add(proteinAcc);
		    }
		   
		    System.out.println(tmpSet);
		    
                    
                    //width=9 height=9 alt="" border="0" hspace="5"><b><a class="genesrch" href="/pathfiles/h_il6Pathway.asp">IL 6 signaling pathway</a></b>
		    
		}
		

	/*	
	    for(int i=1;i<arr.length;i++)
	    {
		System.out.print( arr[i].substring(0, arr[i].indexOf("</a>")) );
		System.out.print("\t");
		System.out.println( "================" );
*/
                
                
	    }
            
	    return;
	}
	catch(Exception e)
	{
	    System.out.println("Error " + e);
	}
    }
}
