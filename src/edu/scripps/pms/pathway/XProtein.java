package edu.scripps.pms.pathway;

import java.util.*;

public class XProtein implements Comparable
{
    private String accession;
    private String id;
    private List<XGene> geneList = new Vector<XGene>();
    private List<String> npList = new Vector<String>();
    private int expression;
    private String[] strArr;

    public XProtein(String accession)
    {
	this.accession = accession;
    }
    
    public XProtein(String accession, int expression)
    {
	this.accession = accession;
	this.expression = expression;
    }

    public XProtein(String accession, int expression, String[] strArr)
    {
	this.accession = accession;
	this.expression = expression;
	this.strArr = strArr;
    }

    public XProtein(String id, String accession, int expression, String[] strArr)
    {
	this(accession, expression, strArr);
	this.id = id;
    }

    public String getId()
    {
	return this.id;
    }
    public XProtein(String accession, XGene gene)
    {
	this.accession = accession;
	addGene(gene);
    }

    public int compareTo(Object obj)
    {
	//    return ((PeptideModel)obj).totalSp - this.totalSp;
	return -1;

    }

    public void addGene(XGene gene)
    {
	this.geneList.add(gene);
    }

    public String getAccession()
    {
	return this.accession;
    }

    public List<XGene> getGeneList()
    {
	return this.geneList;
    }

    public void setExpression(int expression)
    {
	this.expression = expression;
    }

    public int getExpression()
    {
	return this.expression;
    }

    public void addNp(String npName)
    {
	this.npList.add(npName);
    }

    public List<String> getNpList()
    {
	return this.npList;
    }

    public String[] getArr()
    {
	return this.strArr;
    }
    
}
