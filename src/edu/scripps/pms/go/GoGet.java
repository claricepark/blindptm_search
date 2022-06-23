/*
 * @(#)GoGet.java
 *
 * Copyright Notice:
 *
 * Copyright 2004 UBC Bioinformatics Centre
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 *        Free Software Foundation, Inc.,
 *        59 Temple Place, Suite 330, Boston, MA
 *        02111-1307  USA
 * or visit http://www.gnu.org/copyleft/gpl.html
 *
 */

/**
 * @file GoGet.java
 * This is the source file for ubic.atlas.locuslink.db.GoGet
 *
 * @author Tao Xu
 * @date $Date
 */


package edu.scripps.pms.go;

import ubic.db.*;
import ubic.shared.constants.UbicConstants;

import java.sql.*;
import java.util.*;
import java.io.*;





/**
 * @class GoGet
 * This class is a subclass of MysqlDb. It provides an interface to
 * access GO database  
 *
 */
public class GoGet extends MysqlDb {

    public static final String CELLULAR_COMPONENT = "GO:0005575";
    public static final String BIOLOGICAL_PROCESS = "GO:0008150";
    public static final String MOLECULAR_FUNCTION = "GO:0003674";
    public int CELLULAR_COMPONENTTERMID;
    public int BIOLOGICAL_PROCESSTERMID;
    public int MOLECULAR_FUNCTIONTERMID;

    /**
     * @fn GoGet(String dbName)
     * Constructing an GoGet object
     * @param dbName - the name of the database to be accessed
     * @exception throws SQLException
     */
    public GoGet(String dbName) throws SQLException {
        super(dbName);
        getDefulatValues();

    }
    

    /**
     * @fn void closeConnection()
     * Close the DB connection
     * @return void
     * @exception throws SQLException
     */
    public void closeConnection() throws SQLException {
        super.closeConnection();
    }

    public List<String> proteinId2DeepestCellularComponentGoAccs(String pid) throws SQLException {
        return geneProductId2DeepestGoAcc(proteinId2GeneProductId(pid), CELLULAR_COMPONENTTERMID);
    }
    public List<String> proteinId2DeepestMolecularFunctionGoAccs(String pid) throws SQLException {
        return geneProductId2DeepestGoAcc(proteinId2GeneProductId(pid), MOLECULAR_FUNCTIONTERMID);
    }
    public List<String> proteinId2DeepestBiologicalProcessGoAccs(String pid) throws SQLException {
        return geneProductId2DeepestGoAcc(proteinId2GeneProductId(pid), BIOLOGICAL_PROCESSTERMID); 
    }
    public int proteinId2GeneProductId(String pid) throws SQLException {
        int id = -1;
        String q = "SELECT gene_product_id as id FROM  gene_product_synonym WHERE product_synonym = '" + pid + "'";
         
        ResultSet rs = executeQuery(q); 
        if(rs.next()) {
            id = rs.getInt("id");
        }
  
        return id;
    }
    private List<String> geneProductId2DeepestGoAcc(int geneProductId, int termIdOfCategory) 
                                    throws SQLException {
        ArrayList<String> result = new ArrayList<String>();
        String q = "SELECT acc, distance " + 
                   "FROM association a, graph_path g, term t " +
                   "WHERE gene_product_id = " + geneProductId + " AND a.term_id = g.term2_id AND " +
                   "term1_id = " + termIdOfCategory + " AND t.id = a.term_id order by distance DESC";
                  
//System.out.println(q);
        ResultSet rs = executeQuery(q); 
        int maxDistance = -1; 
        while (rs.next()) {
            int distance = rs.getInt("distance");
            if(maxDistance == -1) {
                maxDistance = distance;
            }
            if(distance == maxDistance) {
                result.add(rs.getString("acc"));
            } else {
                break;
            }
        }
  
        return result;
        
    }
    public List<String> prefix2ProteinIds(String prefix) 
                                    throws SQLException {
        ArrayList<String> result = new ArrayList<String>();
        String q = "SELECT product_synonym FROM gene_product_synonym " + 
                   "WHERE product_synonym like '" + prefix + "%'";
                  
//System.out.println(q);
        ResultSet rs = executeQuery(q); 
        while (rs.next()) {
            result.add(rs.getString("product_synonym"));
        }
  
        return result;
        
    }
    private void getDefulatValues() throws SQLException {
        MOLECULAR_FUNCTIONTERMID = goAcc2TermId(MOLECULAR_FUNCTION);
        CELLULAR_COMPONENTTERMID = goAcc2TermId(CELLULAR_COMPONENT);
        BIOLOGICAL_PROCESSTERMID = goAcc2TermId(BIOLOGICAL_PROCESS);
    }
    protected String termId2GoAcc(int termId) throws SQLException {
        String acc = null;
        String q = "SELECT acc FROM  term WHERE id = " + termId;
         
        ResultSet rs = executeQuery(q); 
        if(rs.next()) {
            acc = rs.getString("acc");
        }
        return acc;
    }
    protected int goAcc2TermId(String goAcc) throws SQLException {
        int id = -1;
        String q = "SELECT id FROM  term WHERE acc = '" + goAcc + "'";
         
        ResultSet rs = executeQuery(q); 
        if(rs.next()) {
            id = rs.getInt("id");
        }
        return id;
    }
    public List<String> keyWords2ProteinInfo(String key) throws SQLException {
        LinkedList<String> results = new LinkedList<String>();

        return results;
    }
   
    public int accession2GeneProductId(String ac) throws SQLException {
        int id = -1;
        String q = "SELECT gp.id FROM gene_product gp, gene_product_synonym gps " +
                       "WHERE pgs.gene_product _id = gp.id AND product_synonym = '" + ac + "'";
        ResultSet rs = executeQuery(q); 
        if(rs.next()) {
            id = rs.getInt("gp.id");
        }
        if(id == -1) {

        }
        return id;
    }
    public boolean isAnnotated(String id) throws SQLException {
        String q = "SELECT id FROM dbxref  WHERE xref_key = '" + id + "'";
        ResultSet rs = executeQuery(q); 
        if(rs.next()) {
            return true;
        }
        q = "SELECT gene_product_id as id FROM  gene_product_synonym WHERE product_synonym = '" + id + "'";
        rs = executeQuery(q); 
        if(rs.next()) {
            return true;
        }
        return false;
    }    
    public boolean isSeqInDatabase(String id) throws SQLException {
        String q = "SELECT id FROM seq  WHERE display_id = '" + id + "'";
        ResultSet rs = executeQuery(q); 
        if(rs.next()) {
            return true;
        }
        return false;
    }    
    public void outputAllProteinSeqInFasta(String fileName) throws IOException, SQLException {
        
        PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
        String q = "select seq, description, display_id from seq";
        ResultSet rs = executeQuery(q); 
        while(rs.next()) {
            String seq = rs.getString("seq");
            if(!isNucleicAcidSequence(seq)) {
                // cannot use tab after the identifier, has to be a space
                outFile.println(">" + rs.getString("display_id") + " " + rs.getString("description")); 
                outFile.println(seq); 
            } else {
                System.out.println("nucleic acid found: " + seq);
            }
        }
        
        outFile.close();
    }
    public static boolean isNucleicAcidSequence(String seq) {
        if(seq == null || seq.length() == 0) {
            return false;
        }
        int length = seq.length();
        seq = seq.toUpperCase(); 
        int [] freq = new int[256];
        for(int i = 0; i < length; i++) {
            freq[seq.charAt(i)]++;
        }
          
        return (freq['A'] + freq['T'] + freq['G'] + freq['C'])/length > 0.9; 
        
    }
    public static void main(String [] args) throws Exception {
        //String pid = args[0];
        
        GoGet gg = new GoGet(args[0]);
        String file = args[1];
        String outF = file.substring(0, file.indexOf(".")+1) + "out";
        PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outF)));
        int numProteins = 0;
        int numProteinWithAcc = 0;
        int numAnnotated = 0;
        int numSeqContained = 0; 
        int numPQ = 0;
        BufferedReader br = new BufferedReader(new FileReader(new File(file)));
        String line = br.readLine(); // remove the header line
        String HAHAHA = "hahahaha";
        while((line = br.readLine()) != null) { 
//System.out.println(line);
            String content [] = line.split("\t"); 
            if(content.length < 3) {
                continue;
            }
            
            numProteins++;
            String acc1 = content[2].trim();
          
            String acc2 = HAHAHA;
           
            if(content.length > 6) {
                acc2 = content[7].trim();
            }
            if(acc1.equals("")) {
                acc1 = HAHAHA;
            }
            if(acc2.equals("")) {
                acc2 = HAHAHA;
            }
            if(acc1.equals(HAHAHA) && acc2.equals(HAHAHA)) {
                continue;
            }
            numProteinWithAcc++;
            if(acc1.startsWith("P") || acc1.startsWith("Q") || acc2.startsWith("P") || acc2.startsWith("Q")) {
                numPQ++;
            }
            boolean acc1IsAnnotated = gg.isAnnotated(acc1);
            boolean acc2IsAnnotated =  gg.isAnnotated(acc2);
            boolean isAnnotated = acc1IsAnnotated || acc2IsAnnotated; 
            if(isAnnotated) {
                if(acc1IsAnnotated) {
                    outFile.println(acc1);
                } else {
                    System.out.println(acc2);
                    outFile.println(acc2);
                }
                numAnnotated++;
            }else {
                if(!acc1.equals(HAHAHA)) {
                    outFile.println(acc1); 
                } else {
                    outFile.println(acc2);
                }
     
            }
            boolean isContained = gg.isSeqInDatabase(acc1) || gg.isSeqInDatabase(acc2);
            if(isContained) {
                numSeqContained++;
            }
            System.out.println(acc1 + "\t" + acc2 + "\t is annotated? " + isAnnotated + "\tisSeqIndatabase? " + isContained);
            
            
        }
        System.out.println("NumSeq: " + numProteins +"\tnumProteinWithAcc: " + numProteinWithAcc + "\tnumAnnotated: " + numAnnotated + "\tnumSeqContained: " + numSeqContained + "\tnumPQ: " + numPQ);
        outFile.close();
        /*
        for(String pid : gg.prefix2ProteinIds("IPI")) {
            System.out.println("Deepest GO Annotations for " + pid);
            for(String acc : gg.proteinId2DeepestCellularComponentGoAccs(pid)) {
                System.out.println("Cellular Component: " + acc);
            }        
         
            for(String acc : gg.proteinId2DeepestMolecularFunctionGoAccs(pid)) {
                System.out.println("Molecular Function: " + acc);
            }        
            for(String acc : gg.proteinId2DeepestBiologicalProcessGoAccs("pid")) {
   
                System.out.println("Biological Process: " + acc);
            }        
        }
        */         
    }
}





