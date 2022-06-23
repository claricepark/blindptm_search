/*
* @(#)MysqlDb.java
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
 * @file MysqlDb.java
 * This is the source file for ubic.db.MysqlDb class
 *
 * Author: Tao Xu
 * @date $Date: 2005/05/24 22:52:51 $
 */


package ubic.db;
import java.sql.*;
import ubic.shared.constants.UbicConstants;

/**
 * @class MysqlDb
 * MysqlDb is a subclass of DbConnection. It provides an interface
 * to access MySQL databases.
 *
 * Author: Tao Xu
 * Date: October, 2003
 *
 */
public abstract class MysqlDb extends DbConnection {

    public static final String JDBCDRIVERNAME = "com.mysql.jdbc.Driver";
    public static final String DBSERVERNAME = "mysql";

    /**
    * Constructor that makes the database connection
    * @param dbname the name of the database that has the progdb schema
    * @note This constuctor can be used when the .dbrc file is located at
    * the user's home directory
    */
    protected MysqlDb (String dbname) {

        dbName = dbname;
        jdbcDriverName = JDBCDRIVERNAME;
        dbServerName = DBSERVERNAME;

        makeConnection();

    }



    /**
    * Constructor that makes the database connection
    * @param dbname the name of the database that has the progdb schema
    * @param dbrcFilePath - the path for dbrc file
    */
    protected MysqlDb (String dbname, String dbrcFilePath) {
        dbName = dbname;
        jdbcDriverName = JDBCDRIVERNAME;
        dbServerName = DBSERVERNAME;
        makeConnection(dbrcFilePath);
    }

    /**
     * @fn int getAutoIncrementKey(PreparedStatement prepStmt)
     * Return the auto impremented key of the given PreparedStatement
     * UbicConstants.INVALIDINTID is returned if no auto key is generated
     * @param prepStmt - a PreparedStatement
     * @return An auto incremented key
     * @note This function should be called only after the executeUpdate()
     * function has been called
     * @exception SQLException
     */
    protected int getAutoIncrementKey(PreparedStatement prepStmt)
                                            throws SQLException {

        int id = UbicConstants.INVALIDINTID;

        // get the auto-incremented key
        ResultSet rs = prepStmt.getGeneratedKeys();
        while(rs.next()) {
            id = rs.getInt(1);
        }

        return id;
    }

}
