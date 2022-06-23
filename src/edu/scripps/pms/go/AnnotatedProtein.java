/*
 * @(#)AnnotatedProtein.java
 *
 * Copyright Notice:
 *
 * Copyright 2009 Scripps Research Institute
 *
 *
 */

/**
 * @file AnnotatedProtein.java
 * This is the source file for edu.scripps.pms.go.AnnotatedProtein
 *
 * @author Tao Xu
 * @date $Date
 */


package edu.scripps.pms.go;


import java.util.*;





public class AnnotatedProtein {
    private HashSet<GoTerm> cellularComponent = new HashSet<GoTerm>(5); 
    private HashSet<GoTerm> molecularFunction = new HashSet<GoTerm>(5); 
    private HashSet<GoTerm> biologicalProcess = new HashSet<GoTerm>(5); 
    private String accession;

    public AnnotatedProtein(String acc) {
        accession = acc;
    }

    public HashSet<GoTerm> getMolucularFunctions() {
        return molecularFunction;
    }
    public HashSet<GoTerm> getCellularComponents() {
        return cellularComponent;
    }
    public HashSet<GoTerm> getBiologicalProcesses() {
        return biologicalProcess;
    }

    public void addAnnotation(GoTerm g) {
//System.out.println("GoTerm: " + g);
        if(g.isCellularComponent()) addCellularComponent(g);
        else if(g.isMolecularFunction()) addMolecularFunction(g);
        else if(g.isBiologicalProcess()) addBiologicalProcess(g);
    }
    public void addCellularComponent(GoTerm g) {
        cellularComponent.add(g);
    }
    public void addMolecularFunction(GoTerm g) {
        molecularFunction.add(g);
    }
    public void addBiologicalProcess(GoTerm g) {
        biologicalProcess.add(g);
    }
    
    public String output(Set<String> interestedgoids) {
        StringBuffer sb = new StringBuffer(200);
        sb.append(accession);
        sb.append("\t");
        for(Iterator<GoTerm> it = cellularComponent.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            if(interestedgoids.contains(gt.getGoId())) {
                sb.append(gt.getGoName()); 
                sb.append(';'); 
            }
        }
        sb.append("\t");
        for(Iterator<GoTerm> it = molecularFunction.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            if(interestedgoids.contains(gt.getGoId())) {
                sb.append(gt.getGoName()); 
                sb.append(';'); 
            }
        }
        sb.append("\t");
        for(Iterator<GoTerm> it = biologicalProcess.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            if(interestedgoids.contains(gt.getGoId())) {
                sb.append(gt.getGoName()); 
                sb.append(';'); 
            }
        }

        return sb.toString();
    }
    public String output() {
        StringBuffer sb = new StringBuffer(200);
        sb.append(accession);
        sb.append("\t");
        for(Iterator<GoTerm> it = cellularComponent.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            sb.append(gt.getGoName()); 
            sb.append(';'); 
        }
        sb.append("\t");
        for(Iterator<GoTerm> it = molecularFunction.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            sb.append(gt.getGoName()); 
            sb.append(';'); 
        }
        sb.append("\t");
        for(Iterator<GoTerm> it = biologicalProcess.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            sb.append(gt.getGoName()); 
            sb.append(';'); 
        }

        return sb.toString();
    }
    
}





