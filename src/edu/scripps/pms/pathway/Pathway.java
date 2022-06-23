/*
 * Pathway.java
 *
 * Created on December 27, 2005, 4:25 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.pathway;

import java.util.Set;

/**
 *
 * @author rpark
 */
public class Pathway {
    
    private String name;
    private String url;
    private Set<XGene> geneSet;

    /** Creates a new instance of Pathway */
    public Pathway(String name, String url) {
        this.name = name;
        this.url = url;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
    }
  
    public Set<XGene> getGene()
    {
	return this.geneSet;
    }

    public void setGeneSet(Set geneSet) {
	this.geneSet = geneSet;
    }

    public String toString()
    {
	return this.name;
    }
}
