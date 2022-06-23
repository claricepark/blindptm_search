package edu.scripps.pms.pathway;

import java.util.*;


public class XGene
{
    private String name;
    private List<Pathway> pathwayList = new Vector<Pathway>();
    private List<XProtein> proteinList = new Vector<XProtein>();
    
    public XGene(String name)
    {
	this.setName(name);        
    }

    public void addPathway(String pathwayName, String url)
    {
	pathwayList.add(new Pathway(pathwayName, url));
    }

    public void addProtein(String protein)
    {
	proteinList.add(new XProtein(protein));
    }

    public List<XProtein> getProteinList()
    {
	return this.proteinList;
    }

    public List<Pathway> getPathwayList()
    {
	return this.pathwayList;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString()
    {
	return this.name;
    }
}
