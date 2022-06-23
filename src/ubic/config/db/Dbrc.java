/*
 * @(#)Dbrc.java
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

package ubic.config.db;

import java.io.*;
import java.util.*;
import ubic.config.db.DbrcLoadException;

/**
 * Dbrc class parses .dbrc files stored in the users' home directory.
 * 
 * Author: David He
 * Date:   Jan, 2003
 */
public class Dbrc {

    private static String DBRC_FILE_NAME = ".dbrc";

    private String host;
    private String db;
    private String user;
    private String pwd;

    /**
     * Simple constructor
     */
    public Dbrc() {
	host = "";
	db = "";
	user = "";
	pwd = "";
    }

    /**
     * Constructor that parses the .dbrc file looking for a given
     * database
     * @param dbase the name of the databases that the parser looks
     *              for in the .dbrc file
     */
    public Dbrc(String dbase) throws DbrcLoadException {
	try {
	    load(dbase);
	} catch (DbrcLoadException e) {
	    host = "";
	    db = "";
	    user = "";
	    pwd = "";
	    throw new DbrcLoadException(e.getMessage(), e.dbase, e.dbrcFilePath);
	}
    }

    public Dbrc(String dbase, String dbrcFilePath) throws DbrcLoadException {
	try {
	    loadHelper(dbase, dbrcFilePath);
	} catch (DbrcLoadException e) {
	    throw new DbrcLoadException(e.getMessage(), e.dbase, e.dbrcFilePath);
	}
    }

    public Dbrc(String dbase, String host, String user, String pwd) {
	this.db = dbase;
	this.host = host;
	this.user = user;
	this.pwd = pwd;
    }

    /**
     * This function looks for the .dbrc file in the HOME directory 
     * and if found, the function calls the loadHelper function to 
     * parse the file
     * @param dbase the name of the databases that the parser looks 
     *              for in the .dbrc file
     */
    public void load(String dbase) throws DbrcLoadException {
	String homePath = "";
	String dbrcFilePath = "";

	//get the "HOME" environoment variable 
	homePath = System.getProperty("user.home");

	dbrcFilePath = homePath + "/" + DBRC_FILE_NAME;
	try {
	    loadHelper(dbase, dbrcFilePath);
	} catch (DbrcLoadException e) {
	    throw new DbrcLoadException(e.getMessage(), e.dbase, e.dbrcFilePath);
	}
    }

    /** 
     * This function is similar to the previous load function except that 
     * it takes a dbrcFilePath that points to the .dbrc file instead of looking 
     * at the user's HOME directory by default
     * @param dbase the name of the databases that the parser looks 
     *              for in the .dbrc file
     * @param dbrcFilePath the path to the .dbrc file
     */
    public void load(String dbase, String dbrcFilePath) throws DbrcLoadException {
	try {
	    loadHelper(dbase, dbrcFilePath);
	} catch (DbrcLoadException e) {
	    throw new DbrcLoadException(e.getMessage(), e.dbase, e.dbrcFilePath);
	}
    }

    /**
     * Helper function to retrieve envrionment variables from the system.java
     * This function returns a System.Property object that has the all the 
     * environments variable key value pairs
     */
    private Properties getEnvironment() {

	Properties env = new Properties();
	try {
	    env.load(Runtime.getRuntime().exec("env").getInputStream());
	} catch (IOException e) {
	    System.out.println (e.getMessage());
	}

	return env;
    }

    /**
     * Helper function called by the two load functions. It parses the 
     * dbrc file and looks for the line containing the given database name. When 
     * it finds this line, it sets the db, user, pwd, and host variables of this 
     * class
     * @param dbase the name of the databases that the parser looks 
     *              for in the .dbrc file
     * @param dbrcFilePath the path to the .dbrc file
     */
    private void loadHelper(String dbase, String dbrcFilePath) throws DbrcLoadException {
	String dbrcLine = "";	

	try {
	    FileReader input = new FileReader(dbrcFilePath);
	    BufferedReader buffReader = new BufferedReader(input);

	    while ((dbrcLine = buffReader.readLine()) != null) {
		StringTokenizer tokens = new StringTokenizer(dbrcLine);
		
		if (tokens.hasMoreTokens()) {
		    if (tokens.nextToken().equals(dbase)) {
			this.db = dbase;
			this.user = tokens.nextToken();
			this.pwd = tokens.nextToken();
			this.host = tokens.nextToken();
		    }
		}
	    }

	    input.close();
	} catch (java.io.IOException e) {
	    throw new DbrcLoadException("Can't open dbrc file", dbase, dbrcFilePath);
	}

	if ((this.db == null) || (!this.db.equals(dbase))) {         
	    throw new DbrcLoadException("Can't find database entry", dbase, dbrcFilePath);
	}
    }

    /** 
     * Retrieves the host variable
     */
    public String host() {
	return host;
    }

    /** 
     * Sets the host variable
     * @param host String containing the host variable
     */
    public void host(String host) {
	this.host = host;
    }

    /**
     * Retrieves the user variable
     */
    public String user() {
	return user;
    }

    /** 
     * Sets the user variable
     * @param user String containing the user variable
     */
    public void user(String user) {
	this.user = user;
    }

    /**
     * Retrieves the db variable
     */
    public String db() {
	return db;
    }

    /** 
     * Sets the db variable
     * @param db String containing the db variable
     */
    public void db(String db) {
	this.db = db;
    }

    /**
     * Retrieves the pwd variable
     */
    public String pwd() {
	return pwd;
    }

    /** 
     * Sets the pwd variable
     * @param pwd String containing the pwd variable
     */
    public void pwd(String pwd) {
	this.pwd = pwd;
    }

    /**
     * This function returns a string containing the host, user, pwd and db variables
     */
    public String toString() {
	return (host + " " + user + " " + pwd + " " + db);
    }

    /**
     * This function prints the host, user, pwd and db variables
     */
    public void print() {
	System.out.println (toString());
    }

}



