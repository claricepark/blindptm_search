/*
* @(#)DbConnection.java
*
* Copyright Notice:
*
* Copyright 2001 UBC Bioinformatics Centre
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
 * @file DbConnection.java
 * This is the source file for ubic.db.DbConnection
 *
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:52:51 $
 */


package ubic.db;

import java.sql.*;
import java.util.StringTokenizer;
import ubic.config.db.*;


/**
 * @class DbConnection
 * DBConnection provides connections to databases
 *
 */
public abstract class DbConnection {

    // Connection of this DbConnection
    protected Connection db;
    // Dbrc of this DbConnection
    protected Dbrc dbrc;
    // the name of the database associated with this DbConnection
    protected String dbName;
    // the JDBC Driver name
    protected String jdbcDriverName;
    // the name of the database server
    protected String dbServerName;


    /**
     * @fn void makeConnection(String dbrcFilePath)
     * Make the database connection
     * @param dbrcFilePath - a database name
     * @return void
     */
    protected void makeConnection(String dbrcFilePath) {

        try {
            dbrc = new Dbrc(dbName, dbrcFilePath);
        } catch (DbrcLoadException e) {
            String errMessage = "Error: Can't read .dbrc file located at " +
                    e.dbrcFilePath +
                    " or " +
                    dbName +
                                " entry not found in .dbrc file";

            System.err.println(errMessage);
            System.exit(1);
        }

        // load the Postgresql driver
        try {
            Class.forName(jdbcDriverName);
        } catch (ClassNotFoundException e) {
            System.err.println("No " + dbServerName + " JDBC driver found. " + 
                                           "Check you classpath and try again");
            System.exit(1);
        }

        // Make the connection
        try {

            String connectionURL = "jdbc:" + dbServerName + "://" + 
                                        dbrc.host() + "/" + dbrc.db();

            db = DriverManager.getConnection(connectionURL, dbrc.user(), dbrc.pwd());
        } catch (SQLException e) {
            String errorMessage =
            "Received SQLException while connecting. " + 
            "Check connection parameters and try again" +
            "SQLException: " + e.getMessage();
            System.err.println(errorMessage);
                System.exit(1);
        }
    }


    /**
     * @fn void makeConnection(String host, String user, String pwd)
     * Make the database connection
     * @param host - the name of the database server
     * @param user - user name
     * @param pwd - a password for database
     * @return void
     */
    protected void makeConnection(String host, String user, String pwd) {

        dbrc = new Dbrc(dbName, host, user, pwd);
        // load the Postgresql driver
        try {
            Class.forName(jdbcDriverName);
        } catch (ClassNotFoundException e) {
            System.err.println("No " + dbServerName + " JDBC driver found. " + 
                                            "Check you classpath and try again");
            System.exit(1);
        }

        // Make the connection
        try {

            String connectionURL = "jdbc:" + dbServerName + "://" + dbrc.host() + "/" + dbrc.db();


            db = DriverManager.getConnection(connectionURL, dbrc.user(), dbrc.pwd());
        } catch (SQLException e) {
            String errorMessage =
            "Received SQLException while connecting. Check connection parameters and try again" +
            "SQLException: " + e.getMessage();
            System.err.println(errorMessage);
            System.exit(1);
        }
    }

    /**
     * @fn String getDbName()
     * Return the database name of this DBConnection object
     * @return A String for database name
     */
    public String getDbName() {

        return dbName;
    }



    /**
     * @fn boolean isClosed()
     * Return true if the Connection object db of this DbConnection is closed,
     * otherwise return false
     * @return true if this DbConnection is closed.
     * @exception throws SQLException
     */
    public boolean isClosed() throws SQLException {

        return db.isClosed();

    }

    /**
     * @fn void makeConnection()
     * Makes a connection to the database
     * @return void
     */
    protected void makeConnection() {

        // get database information from dbrc file
        try {
            dbrc = new Dbrc(dbName);
        } catch (DbrcLoadException e) {

            String errMessage = "Error: Can't read .dbrc file located at " +
                                e.dbrcFilePath + " or " + dbName +
                                " entry not found in .dbrc file";

            System.err.println(errMessage);
            System.exit(1);
        }

        // load the Postgresql driver
        try {
            Class.forName(jdbcDriverName);
        } catch (ClassNotFoundException e) {
            System.err.println("No " + dbServerName +
                " JDBC driver found. Check you classpath and try again");
            System.exit(1);
        }

        // Make the connection
        try {

            String connectionURL = "jdbc:" + dbServerName + "://" +
                            dbrc.host() + "/" + dbrc.db();

            db = DriverManager.getConnection(
                connectionURL, dbrc.user(), dbrc.pwd());

        } catch (SQLException e) {

            String errorMessage =
                "Received SQLException while connecting. " +
                "Check connection parameters and try again" +
                "SQLException: " + e.getMessage();

            System.err.println(errorMessage);
            System.exit(1);
        }
    }

    /**
     * @fn ResultSet executeQuery(String query)
     * Returns a ResultSet object containing the result of running the given
     * mysql query
     * @param query - String containg the mysql query.
     * @return A ResultSet
     * @exception SQLException
     */
    protected ResultSet executeQuery(String query) throws SQLException {

        return db.createStatement().executeQuery(query);
    }

    /**
     * @fn Statement createStatement()
     * Creates a statement to query the current database
     * @return A Statement
     * @exception SQLException
     */
    protected Statement createStatement() throws SQLException {
        return db.createStatement();
    }



    /**
     * @fn PreparedStatement prepareStatement(String query)
     * Returns a PreparedStatement associated with given query string
     * @param query - a mysql query String.
     * @return A PreparedStatement object
     * @exception SQLException
     */
    protected PreparedStatement prepareStatement(String query) throws SQLException {

        return db.prepareStatement(query);
    }

    /**
     * @fn static int getNumRows(ResultSet rs)
     * Returns the number of rows in the given ResultSet object
     * @param rs - A ResultSet
     * @return Number of rows stored in the given ResultSet object.
     * @note the cursor of the ResultSet object is reset to beforeFirst after
     * the execution of this method
     * @exception SQLException
     */
    protected static int getNumRows(ResultSet rs) throws SQLException{

        rs.last();
        int numRows = rs.getRow();
        rs.beforeFirst();

        return numRows;
    }

    /**
     * @fn static boolean isEmptyResultSet(ResultSet rs)
     * Returns true if the number of rows in the given ResultSet object is 0,
     * otherwise return false
     * @param rs - A ResultSet
     * @return - true if the ResultSet is empty.
     * @exception SQLException
     */
    protected static boolean isEmptyResultSet(ResultSet rs) throws SQLException {

        return getNumRows(rs) == 0;
    }


    /**
     * @fn String escapeSqlWildCards(String targetString)
     * Return a string that escapes the SQL wild cards "_" and "%" in
     * the input string stringToBeEscaped
     * @param targetString - the string to be escaped
     * @return the escaped string object
     * @note Only the string appears in LIKE CAN to be escaped,
     * @note DO NOT USE THIS FUNCTION WHEN EXACT MATCH IS DESIRED
     */
    protected String escapeSqlWildCards(String targetString) {

        if (targetString == null) {
            return targetString;
        }

        // escaping "_" character
        String result = targetString.replaceAll("_", "\\\\_");

        // escaping "%" character
        result = result.replaceAll("%", "\\\\%");


        return result;

    }


    /**
     * @fn void closeConnection()
     * Closes the connection the database
     * @return void
     * @exception SQLException
     */
    protected void closeConnection() throws SQLException {
        db.close();
    }

}
