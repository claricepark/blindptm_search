/*
* @(#)PostgresqlDb.java
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
 * @file PostgresqlDb.java
 * This is the source file for ubic.db.PostgresqlDb class
 *
 * Author: Tao Xu
 * @date $Date: 2005/05/24 22:52:51 $
 */

package ubic.db;
import java.sql.*;
import java.io.*;


/**
 * @class PostgresqlDb
 * PostgresqlDb provides connections to PostgreSQL databases
 *
 */
public abstract class PostgresqlDb extends DbConnection {

    // the postgresql jdbc driver name
    public static final String JDBCDRIVERNAME = "org.postgresql.Driver";

    // the postgresql database server name
    public static final String DBSERVERNAME = "postgresql";

    /**
     *@fn PostgresqlDb (String dbname)
     * Constructor that makes the database connection
     * @param dbname - the name of the database
     */
    protected PostgresqlDb (String dbname) {

        dbName = dbname;
        jdbcDriverName = JDBCDRIVERNAME;
        dbServerName = DBSERVERNAME;
        makeConnection();

    }


    /**
     *@fn PostgresqlDb (String dbname, String dbrcFilePath)
     * Constructor that makes the database connection
     * @param dbname - the name of the database
     * @param dbrcFilePath - the path for the .dbrc file
     */
    protected PostgresqlDb (String dbname, String dbrcFilePath) {

        dbName = dbname;
        jdbcDriverName = JDBCDRIVERNAME;
        dbServerName = DBSERVERNAME;
        makeConnection(dbrcFilePath);
    }


    /**
     *@fn PostgresqlDb (String dbname, String dbrcFilePath)
     * Constructor that makes the database connection
     * @param dbname - the name of the database
     * @param host - the host of the database server
     * @param user - user name
     * @param pwd - password
     */
    protected PostgresqlDb (String dbname, String host, String user, String pwd) {
        dbName = dbname;
        jdbcDriverName = JDBCDRIVERNAME;
        dbServerName = DBSERVERNAME;
        makeConnection(host, user, pwd);

    }
}
