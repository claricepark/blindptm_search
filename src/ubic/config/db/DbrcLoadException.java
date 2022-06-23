/*
 * @(#)DbrcLoadException.java
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

/**
 * DbrcLoadException extends the Exception class and is thrown when 
 * the .dbrc file is not found, when it can't be opend, or when the 
 * database name is not found in the .dbrc file
 */
public class DbrcLoadException extends Exception {

    /** Explains the cause of the error */
    public String errorType;
    public String dbase;
    public String dbrcFilePath;

    /** 
     * Constructor that calls the constructor of the super class and sets 
     * the necessary variables
     */
    public DbrcLoadException (String errorType, String dbase, String dbrcFilePath) {
	super (errorType);
	this.dbase = dbase;
	this.dbrcFilePath = dbrcFilePath;
    }
}





