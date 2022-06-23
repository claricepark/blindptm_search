

/**
 * @author Tao Xu 
 * @version $Id: InteractionGet.java,v 1.1 2005/05/09 23:40:32 taoxu Exp $
 * @date $Date: 2005/05/09 23:40:32 $
 */


package edu.scripps.pms.interaction.db; 

import ubic.db.*;
import ubic.shared.constants.UbicConstants;

import java.sql.*;
import java.util.*;



import org.apache.log4j.*;

/**
 * @class InteractionGet
 * This class is a subclass of MysqlDb. It provides an interface to
 * access Atlas HomoloGene database.
 *
 */
public class InteractionGet extends MysqlDb {



    /**
     * @fn InteractionGet(String dbName)
     * Constucting a InteractionGet object
     * @param dbName - the name of the database to be accessed
     * @exception SQLException
     */
    public InteractionGet(String dbName) throws SQLException {
        super(dbName);

    }
   
    private int taxonId2OrganismId(int taxonId) throws SQLException {
        int organismId = -1;
        String query = "SELECT " +
                           "organism_id " +
                       "FROM " +
                           "Organism " +
                       "WHERE " +
                           "taxonid = " + taxonId;
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            organismId = rs.getInt("organism_id");
        }   
        return organismId;
    }
    
    public boolean isIncluded(String shortLabel) throws SQLException {
        boolean result = false;
        String query = "SELECT " +
                           "interactor_id " + 
                       "FROM " +
                           "Interactor " +
                       "WHERE " +
                           "shortlabel = '" + shortLabel + "'";
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            result = true;
            break;
        }   
        return result;
    } 
    public boolean isIncluded(String shortLabel, int taxonId) throws SQLException {
        boolean result = false;
        int organismId = taxonId2OrganismId(taxonId);
        String query = "SELECT " +
                           "interactor_id " + 
                       "FROM " +
                           "Interactor " +
                       "WHERE " +
                           "shortlabel = '" + shortLabel + "' AND " +
                           "organism_id = " + organismId;
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            result = true;
            break;
        }   
        return result;
    } 
    public Set<String> shortLabel2InteractionPartners(String label)
                         throws SQLException {
        Set<String> result = new HashSet<String>();
        String query = "SELECT " +
                           "i2.shortlabel " +
                       "FROM " +
                           "Interactor i1, Interactor i2, Interactor_Interaction ii1, " +
                           "Interactor_Interaction ii2 " +
                       "WHERE " +
                           "i1.shortlabel = '" + label + "' AND " +
                           "ii1.interactor_id = i1.interactor_id AND " +
                           "ii2.interaction_id = ii1.interaction_id AND " +
                           "i2.interactor_id = ii2.interactor_id AND " +
                           "i2.interactor_id <> i1.interactor_id AND " +
                           "i1.organism_id = i2.organism_id"; 
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            result.add(rs.getString("shortlabel").toUpperCase());
        }   

        return result;

    }
    public Set<String> shortLabel2InteractionPartners(String label, int taxonId)
                         throws SQLException {
        int organismId = taxonId2OrganismId(taxonId);
        Set<String> result = new HashSet<String>();
        String query = "SELECT " +
                           "i2.shortlabel " +
                       "FROM " +
                           "Interactor i1, Interactor i2, Interactor_Interaction ii1, " +
                           "Interactor_Interaction ii2 " +
                       "WHERE " +
                           "i1.shortlabel = '" + label + "' AND " +
                           "ii1.interactor_id = i1.interactor_id AND " +
                           "ii2.interaction_id = ii1.interaction_id AND " +
                           "i2.interactor_id = ii2.interactor_id AND " +
                           "i2.interactor_id <> i1.interactor_id AND " +
                           "i1.organism_id = i2.organism_id AND " +
                           "i1.organism_id = " + organismId;
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            result.add(rs.getString("shortlabel").toUpperCase());
        }   

        return result;

    }

    protected List<Integer> shortLabel2InteractorIds(String label) throws SQLException {
        List<Integer> result = new ArrayList<Integer>();

        String query =  "SELECT " +
                            "interactor_id " +
                        "FROM " +
                            "Interactor " +
                        "WHERE " +
                            "shortlabel = '" + label + "'";

        ResultSet rs = executeQuery(query);


        while(rs.next()) {
            result.add(new Integer(rs.getInt("interactor_id")));
        }
        
        return result;
    }
   
    public static void main(String [] args) {
        try {
            String shortLabel = "YPR174C";
            InteractionGet ig = new InteractionGet("bind");
            List<Integer> ids = ig.shortLabel2InteractorIds(shortLabel);
            for (Integer it : ids) {
                System.out.println(it);
            }
            // 4932 is taxonId for yeast 
            for(String label : ig.shortLabel2InteractionPartners(shortLabel, 4932)) {
                System.out.println(label);
            }
            ig.closeConnection();
        } catch(SQLException e) {
            e.printStackTrace();
        }
    }

    /**
     * @fn void closeConnection()
     * Close the DB connection
     * @return void
     * @exception SQLException
     */
    public void closeConnection() throws SQLException {
        super.closeConnection();
    }

}





