package edu.scripps.pms.pathway;

import java.io.BufferedReader;
import java.util.*;
import java.io.*;
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
 * @version $Id: MapFileGenerator.java,v 1.5 2007/03/05 23:09:22 rpark Exp $
 */
public class MapFileGenerator 
{
    public static void main(String args[]) throws IOException, Exception
    {
	MapFileGenerator gen = new MapFileGenerator();
	gen.generateMapFiles(args[0]);
    }

    public void generateMapFiles(String fileName) throws IOException, Exception
    {
	//Generate gene to protein map file
	//Generate protein to gene map file
	//Generate gene to pathway map file
	Hashtable<String, XProtein> ht = MapFileGenerator.getXRef(fileName);
	Hashtable<String, Set> gpTable = new Hashtable<String, Set>();

	FileOutputStream out = new FileOutputStream(fileName + ".pg.map");
	PrintStream p = new PrintStream( out );

	Set geneSet = new HashSet();

	for (Enumeration<String> e = ht.keys(); e.hasMoreElements(); )
	{
	    String key = e.nextElement();
	    XProtein protein = ht.get(key);
	//    System.out.print(protein.getAccession() + "\t");
	    p.print(protein.getAccession() + "\t");
	    List l = new Vector();

	    List<XGene> gList = protein.getGeneList();
	    for(Iterator<XGene> itr=gList.iterator(); itr.hasNext(); )
	    {
		XGene gene = itr.next();
		if( !l.contains(gene.getName()) )
		{
		    l.add(gene.getName());
		    geneSet.add(gene.getName());
		    p.print(gene.getName());
		    p.print(";");
		}

		//generate gene to protein map file
		Set tempSet = gpTable.get(gene.getName());
		if(null == tempSet)
		{
		    tempSet = new HashSet();
		    tempSet.add(protein.getAccession());
		    gpTable.put(gene.getName(), tempSet);
		}
		else
		{
		    tempSet.add(protein.getAccession());
		}
	    }

	    //XRef ref = new XRef(); 
	    //ref.readBioCarda( gList, protein.getAccession() );
	    p.println();
	}

	p.close();
	out.close();

	out = new FileOutputStream(fileName + ".gp.map");
	p = new PrintStream( out );

	for (Enumeration<String> e = gpTable.keys(); e.hasMoreElements(); )
	{
	    String geneName = e.nextElement();
	    Set set = gpTable.get(geneName);

	    p.print(geneName); p.print('\t');

	    for(Iterator<String> itr=set.iterator(); itr.hasNext(); )
	    {
		String pStr = itr.next();
		p.print(pStr); p.print(';');
	    }

	    p.println();
	}

	p.close();
	out.close();


	File f = new File(fileName + ".pathway.map");

	if( f.exists() )
	{

	    String eachLine;
	    BufferedReader br = new BufferedReader(new FileReader(f));

	    Set set = new HashSet();
	    while( null != (eachLine = br.readLine()) )
	    {
		String[] arr = eachLine.split("\t");
		set.add(arr[0]);
	    }

	    System.out.println("appending to existing pathway file..");

	    BufferedWriter bw = new BufferedWriter(new FileWriter(f, true));

	    for(Iterator<String> itr = geneSet.iterator(); itr.hasNext(); )
	    {
		String geneName = itr.next();

		if( set.contains(geneName) )
		    continue;

		bw.write(geneName); bw.write('\t');
		bw.write( XRef.readBioCarda(geneName) );
		bw.newLine();
		bw.flush();

		System.out.print(".");
	    }

	    bw.close();
	}
	else
	{
	    //out = new FileOutputStream(fileName + ".pathway.map");
	    out = new FileOutputStream(f);
	    p = new PrintStream( out );
	    System.out.println("generating pathway");
	    for(Iterator<String> itr = geneSet.iterator(); itr.hasNext(); )
	    {
		String geneName = itr.next();
		p.print(geneName); p.print('\t');
		p.println( XRef.readBioCarda(geneName) );
		System.out.print(".");
	    }

	}



	p.close();
	out.close();

	//generatePathwayToGeneMap(fileName + ".pathway.map");
    }

    public static Hashtable<String, XProtein> getXRef(String fileName) throws IOException, Exception
    {
	Hashtable<String, XProtein> ht = new Hashtable<String, XProtein>();
	String lastLine = "";

        BufferedReader br = new BufferedReader(new FileReader(fileName));

	lastLine = br.readLine();  //read out comment

	while( (lastLine=br.readLine())!= null )
	{
            String[] str = lastLine.split("\t");

	    
	    //System.out.println(str[5] + "=="  + str[6]);
	    //System.out.println(str[7] + "=="  + str[8] + " " + str[9]);


	    

	    //for HUMAN and MOUSE
	    /*
	    String[] pArr = str[7].split(";");

	    for(int i=0;i<pArr.length;i++)
	    {
		//String geneName;
		if( !"".equals(str[6]) )
		    addGene(ht, pArr[i], str[6]);

		if( !"".equals(str[5]) )
		    addGene(ht, pArr[i], str[5]);
	    }	    
	    */

	    //For RAT, they have different column #

	    if(str.length<=9)
		continue;

	    String[] pArr = str[9].split(";");

	    for(int i=0;i<pArr.length;i++)
	    {
		//String geneName;
		if( !"".equals(str[7]) )
		    addGene(ht, pArr[i], str[7]);

		if( !"".equals(str[8]) )
		    addGene(ht, pArr[i], str[8]);
	    }	    

	    if(str.length<=10)
		continue;


	    //for mouse
	    String[] npArr = str[7].split(";");
	    for(int i=0;i<npArr.length;i++)
	    {
		//String geneName;
		if( !"".equals(str[10]) )
		{
		    XProtein xprotein = ht.get(npArr[i]);

		    if(null == xprotein)
		    {
			xprotein = new XProtein(npArr[i]);
			xprotein.addNp(str[10]);
			ht.put(npArr[i], xprotein);
		    }
		    else
			xprotein.addNp(str[10]);
			

		}
	    }	    
	}

	return ht;
    }

    public static void addNp(Hashtable<String, XProtein> ht, String accession, String geneName)
    {

    }
    
    public static void addGene(Hashtable<String, XProtein> ht, String accession, String geneName)
    {
	geneName = geneName.substring( geneName.indexOf(',')+1 );

	XProtein temp = ht.get(accession);

	if(null == temp)
	{
	    ht.put(accession, new XProtein(accession, new XGene(geneName) ));
	}
	else
	{
	    temp.addGene( new XGene(geneName) );
	    //    System.out.println(temp.getAccession() + " " + temp.getGeneList().size());
	}
    }
	
}
