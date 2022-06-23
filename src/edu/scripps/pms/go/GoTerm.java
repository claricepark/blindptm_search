/*
 * @(#)GoTerm.java
 *
 * Copyright Notice:
 *
 * Copyright 2009 Scripps Research Institute
 *
 *
 */

/**
 * @file GoTerm.java
 * This is the source file for edu.scripps.pms.go.GoTerm
 *
 * @author Tao Xu
 * @date $Date
 */


package edu.scripps.pms.go;


import java.util.*;





public class GoTerm {
    private char goType = ' ';   // c for cellular, b for biological, m for molecular function
    private String goId;
    private String goName;

    // goids of parent terms, including both is_a and relationships such as has_part, part_of, regulates, positively_regulates, negatively_regulates
    // also need to consider intersection_of, need to think of has_part
    //private ArrayList<String> isas = new ArrayList<String>(); // goids of parent terms, including both is_a and relationships
    private HashSet<String> isas = new HashSet<String>(); // goids of parent terms, including both is_a and relationships

    public GoTerm(String l) {
        String [] arr = l.split("\t");        
        if(arr != null && arr.length > 3) { 
            goType = arr[2] == null || arr[2].length() < 1? ' ' : arr[2].charAt(0);
            goId = arr[3];
            goName = arr[1];
        } else {System.out.println("Strange term: " + l); }
    }
   
    public GoTerm() {}
 
    public GoTerm(char c, String id, String term) {
        goType = c;
        goId = id;
        goName = term;

    }

    // get parent terms, including is_a and relationships
    public Iterator<String> getIsAs() {
        return isas.iterator();
    }
    public void addIsA(String is_a) {
        isas.add(is_a);
    }

    // store the relations in isas
    public void addRelationship(String relationship) {
        isas.add(relationship);
    }
    // store the intersection_of in isas
    public void addIntersectionOf(String s) {
        isas.add(s);
    }
    public void setGoId(String id) {
        goId = id;
    }
    public String getGoId() {
        return goId;
    }
   
    public void setGoName(String name) {
        goName = name;
    } 
    public String getGoName() {
        return goName;
    }

    public void setGoType(char type) {
        goType = type;
    } 
    public char getGoType() {
        return goType;
    }

    
    public boolean isMolecularFunction() {
        return goType == 'm';
    }
    public boolean isCellularComponent() {
        return goType == 'c';
    }
    public boolean isBiologicalProcess() {
        return goType == 'b';
    }
}





