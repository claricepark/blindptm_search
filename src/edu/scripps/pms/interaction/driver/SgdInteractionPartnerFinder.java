

/**
 * @author Tao Xu 
 * @version $Id
 * @date $Date
 */



import edu.scripps.pms.interaction.db.InteractionGet;
import java.sql.*;
import java.util.*;
import java.io.*;


import org.apache.log4j.*;

public class SgdInteractionPartnerFinder {
    // assume the interactions.tab from SGD is in the current working directory

    private static final String interactionFile = "interactions.tab";
    private static final String USAGE = "Usage: java SgdInteractionPartnerFinder proteinListFile";
    private static final String format = "Locus\tIsInDatabase\tNumPartnersInList\tPartnersInList\tPartnersInDatabase";
    private HashSet<String> proteins = new HashSet<String>();
    private LinkedList<String> proteinList = new LinkedList<String>();
    private String outputFile;
    
    private HashMap<String, HashSet<String>> sgdInteractions = new HashMap<String, HashSet<String>>();
   
    public SgdInteractionPartnerFinder(String proteinListFile)
                throws SQLException, IOException {
        int pos = proteinListFile.indexOf(".");
        outputFile = proteinListFile.substring(0, pos) + ".out";
        loadInteractions();
        loadProteins(proteinListFile);
    }
    public static void main(String [] args) {
        try {
            SgdInteractionPartnerFinder ipf = new SgdInteractionPartnerFinder(args[0]);
            //int taxonId = Integer.parseInt(args[1]);
            ipf.outputResult();
        } catch(Exception e) {
            System.err.println(USAGE);
            e.printStackTrace();
        }
    }
    public void outputResult() throws SQLException, IOException {
        PrintStream ps = new PrintStream(outputFile);
        ps.println(format);
        int numProteins = 0;
        int numProteinsInDatabase = 0;
        int numProteinsWithPartnersInList = 0;
        int numProteinsNotInDatabase = 0;
        
        for (String shortLabel : proteinList) {
            numProteins++;
//System.out.println("shortLabel: " + shortLabel);
            String output = shortLabel;
            if (isIncluded(shortLabel)) {
                numProteinsInDatabase++;
                output += "\tY"; 
                Set<String> partners = id2InteractionPartners(shortLabel);
                String allPartners = getAllPartners(partners);
                partners.retainAll(proteins);
                if(partners.size() > 0) {
                    numProteinsWithPartnersInList++;
                }
                String partnersIncluded = getAllPartners(partners);
                output += "\t" + partners.size() + "\t" + partnersIncluded + "\t" + allPartners;
            } else {
                numProteinsNotInDatabase++;
                output += "\tN";
System.out.println("Protein not in sgdInteraction.tab file: " + shortLabel);
            }
//System.out.println(output);
            ps.println(output);
        }
        ps.println("Number of proteins identified: " + numProteins);
        ps.println("Number of identified proteins included in the database: " + numProteinsInDatabase);
        ps.println("Number of identified proteins not included in the database: " + numProteinsNotInDatabase);
        ps.println("Number of identified proteins with interacting partners found: " + numProteinsWithPartnersInList);
        ps.close();
    }
    public Set<String> id2InteractionPartners(String id) {
        return new HashSet<String>(sgdInteractions.get(id));
    }
    public boolean isIncluded(String id) {
        return sgdInteractions.get(id) != null;
    }
    public String getAllPartners(Set<String> partners) {
        String output = "";
        if (partners.size() == 0) {
            output = "NONE";
        } 
        for (String label : partners) {
            output += label + " ";
        }

        return output.trim();
    }
    private void loadInteractions() throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(interactionFile));
        String line = br.readLine().trim();
       
        while (line != null) {
            String [] elements = line.split("\t");
            String [] genes = elements[1].split("\\|");
            if(genes.length > 2) {
                System.out.println("more than one interactor ids found: " + genes.length);
                System.out.println("Line: " + elements[1]);
            } else {
                String id1 = genes[0].split(" ")[0];
                String id2 = genes[1].split(" ")[0];
//                System.out.println("Line: " + elements[1] + "\t id1: " + id1 + "\tid2: " + id2);
                addInteraction(id1, id2);
            }
            
            line = br.readLine();
        }
        System.err.println("Number of Proteins in the set: " + sgdInteractions.size());
        br.close();
    }
    private void addInteraction(String id1, String id2) {
        HashSet<String> id1Set = sgdInteractions.get(id1);
        HashSet<String> id2Set = sgdInteractions.get(id2);
        if(id1Set != null) {
            id1Set.add(id2);
        } else {
            id1Set = new HashSet<String>();
            id1Set.add(id2);
            sgdInteractions.put(id1, id1Set);
        }
        if(id2Set != null) {
            id2Set.add(id1);
        } else {
            id2Set = new HashSet<String>();
            id2Set.add(id1);
            sgdInteractions.put(id2, id2Set);
        }
    }
    private void loadProteins(String proteinListFile) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(proteinListFile));
        String line = br.readLine().trim().toUpperCase();
        while (line != null) {
            proteinList.add(line);
            proteins.add(line);
            line = br.readLine();
        }
        //System.err.println("Number of Proteins in the set: " + proteins.size());
        br.close();
    }
}





