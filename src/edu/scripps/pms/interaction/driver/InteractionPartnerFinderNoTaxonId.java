

/**
 * @author Tao Xu 
 * @version $Id: InteractionPartnerFinderNoTaxonId.java,v 1.1 2005/05/09 23:41:34 taoxu Exp $
 * @date $Date: 2005/05/09 23:41:34 $
 */



import edu.scripps.pms.interaction.db.InteractionGet;
import java.sql.*;
import java.util.*;
import java.io.*;


import org.apache.log4j.*;

public class InteractionPartnerFinderNoTaxonId {

    private static final String USAGE = "Usage: java InteractionPartnerFinderNoTaxonId proteinListFile taxonId";
    private static final String format = "Locus\tIsInDatabase\tNumPartnersInList\tPartnersInList\tPartnersInDatabase";
    private HashSet<String> proteins = new HashSet<String>();
    private LinkedList<String> proteinList = new LinkedList<String>();
    private InteractionGet ig;
    private String outputFile;
   
    public InteractionPartnerFinderNoTaxonId(String dbName, String proteinListFile)
                throws SQLException, IOException {
        ig = new InteractionGet(dbName);
        int pos = proteinListFile.indexOf(".");
        outputFile = proteinListFile.substring(0, pos) + ".out";
        loadProteins(proteinListFile);
    }
    public static void main(String [] args) {
        try {
            InteractionPartnerFinderNoTaxonId ipf = new InteractionPartnerFinderNoTaxonId("bind", args[0]);
            //int taxonId = Integer.parseInt(args[1]);
            ipf.outputResult();
            ipf.closeDbConnection();
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
            String output = shortLabel;
            if (ig.isIncluded(shortLabel)) {
                numProteinsInDatabase++;
                output += "\tY"; 
                Set<String> partners = ig.shortLabel2InteractionPartners
                                (shortLabel);
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
            }
            ps.println(output);
        }
        ps.println("Number of proteins identified: " + numProteins);
        ps.println("Number of identified proteins included in the database: " + numProteinsInDatabase);
        ps.println("Number of identified proteins not included in the database: " + numProteinsNotInDatabase);
        ps.println("Number of identified proteins with interacting partners found: " + numProteinsWithPartnersInList);
        ps.close();
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
    public void closeDbConnection() throws SQLException {
        ig.closeConnection();
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





