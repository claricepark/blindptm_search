/*
 * @(#)GeneOntologyExtReader.java
 *
 * Copyright Notice:
 *
 * Copyright 2010 Scripps Research Institute
 *
 *
 */

/**
 * @file GeneOntologyExtReader.java
 * This is the source file for edu.scripps.pms.go.GeneOntologyExtReader
 *
 * @author Tao Xu
 * @date $Date
 */


package edu.scripps.pms.go;

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.HashSet; 

public class GeneOntologyExtReader
{
    private ArrayList<GoTerm> goTerms=new ArrayList<GoTerm>();
    private String goFile;

    public GeneOntologyExtReader(String filename) throws IOException
    {
        goFile=filename;
    }

    // read GoTerm from Go.defs file
    public ArrayList<GoTerm> getGoTerms() throws IOException
    {
        BufferedReader br= new BufferedReader(new FileReader(goFile));
        String line= br.readLine();
        StringBuffer termLine= new StringBuffer();
        
        while(line!=null) //after goid: GO:0050369
        {
                if(line.startsWith("term:"))
                        termLine.append("\t"+line.substring(6));
                else if(line.startsWith("goid:"))
                        termLine.append("\t\t"+line.substring(6));
                else if(line.equals(""))
                {
                        goTerms.add(new GoTerm(termLine.toString()));//how to check if repeated?
                        termLine=new StringBuffer();
                }
                line=br.readLine();
        }
        br.close();			

        return goTerms;
    }


    // read GoTerm from gene_ontology_ext.obo file
    public ArrayList<GoTerm> getExtGoTerms() throws IOException
    {
        BufferedReader br= new BufferedReader(new FileReader(goFile));
        String line= br.readLine();
        while(line != null && !line.startsWith("[Term")) line = br.readLine(); // read till the firt term 
        GoTerm goterm = new GoTerm();
                        

        while(line!=null) //after goid: GO:0050369
        {
//System.out.println(line);
            if(line.startsWith("[Term")) {
                goterm = new GoTerm();
            } else if(line.startsWith("id:")) {
                goterm.setGoId(line.substring(4)); 
            } else if(line.startsWith("name:")) {
                goterm.setGoName(line.substring(6)); 
            } else if(line.startsWith("namespace:")) {
                goterm.setGoType(line.charAt(11)); 
            } else if(line.startsWith("is_a:")) {
                goterm.addIsA(line.split(" ")[1]); 
            } else if(line.startsWith("relationship:")) { // also need to consider intersection_of
                goterm.addRelationship(line.split(" ")[2]); 
            } else if(line.startsWith("intersection_of:")) { // also need to consider intersection_of
                String [] temp1 = line.split(" ");
                if(temp1[1] != null && temp1[1].startsWith("GO")) {

                    //intersection_of: GO:0065007 ! biological regulation
                    goterm.addRelationship(temp1[1]); 
                } else if(temp1.length > 3) {
                    //intersection_of: regulates GO:0009058 ! biosynthetic process
                    if(temp1[2] != null && temp1[2].startsWith("GO")) {
                        goterm.addRelationship(temp1[2]); 
                    }
                }
            } else if(line.equals("")) {
                goTerms.add(goterm);//how to check if repeated?
            } else if(line.startsWith("[Typedef")) { // ignore Typedefs
                break;
            }
            line=br.readLine();
        }
        br.close();			

        return goTerms;
    }
    
    public static void main(String [] args) throws IOException {
        GeneOntologyExtReader goer = new GeneOntologyExtReader("/home/taoxu/taoxu_on_data/projects/goa/20100326/gene_ontology_ext.obo");
        ArrayList<GoTerm> terms = goer.getExtGoTerms();
        Iterator<GoTerm> it = terms.iterator();
        while(it.hasNext()) {
            GoTerm goterm = it.next();
            StringBuffer sb = new StringBuffer();
            Iterator<String> pit = goterm.getIsAs();
            while(pit.hasNext()) sb.append(pit.next() +"|");
            System.out.println("goid:" + goterm.getGoId() + "\tgoname:" + goterm.getGoName() + "\tgotype:" + goterm.getGoType() + "\tis_as:" + sb);
        }
    }
        
}
