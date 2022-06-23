/*
 * @(#)Dbxref.java
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
 * @file Dbxref.java
 * Class file for Dbxref class
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.io.*;
import java.util.*;
import org.apache.log4j.*;


public class Dbxref {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Dbxref.class);

    // Private members of Dbxref class
    private int dbxrefId = UbicConstants.INVALIDINTID;
    private String idType = null;
    private String idValue = null;
    private String dbSource = null;



    /**
     * @fn Dbxref (int dbxrefId, String idType, String idValue, String dbSource)
     * Constructor for Interactor
     * @param dbxrefId - a dbxref_id
     * @param idType - idtype
     * @param idValue - idvalue
     * @param dbSource - dbsource
     */
    public Dbxref (int dbxrefId, String idType, String idValue, String dbSource){
        this.dbxrefId = dbxrefId;
        this.idType = idType;
        this.idValue = idValue;
        this.dbSource = dbSource;
    }


    /**
     * @fn Dbxref (int String idType, String idValue, String dbSource)
     * Constructor for Interactor
     * @param idType - idtype
     * @param idValue - idvalue
     * @param dbSource - dbsource
     */
    public Dbxref (String idType, String idValue, String dbSource) {

        this.idType = idType;
        this.idValue = idValue;
        this.dbSource = dbSource;
    }

    // Getters

    /**
     * @fn int getDbxrefId()
     * Return the dbxref id of this Dbxref.
     * @return the dbxref id.
     */
    public int getDbxrefId() {
        return dbxrefId;
    }

    /**
     * @fn String getIdType()
     * Return the id type of this Dbxref.
     * @return the id type.
     */
    public String getIdType() {
        return idType;
    }



    /**
     * @fn String getIdValue()
     * Return the id value of this Dbxref.
     * @return the id value.
     */
    public String getIdValue() {
        return idValue;
    }



    /**
     * @fn String getDbSource()
     * Return the db source of this Dbxref.
     * @return the db source.
     */
    public String getDbSource() {
        return dbSource;
    }

    /**
     * @fn String toString() 
     * Return a string that represent this Dbxref object.
     * @return a string that represent this Dbxref object.
     */
    public String toString() {
        return "dbSource: " + dbSource + "\tidType: " + 
                        idType + "\tidValue: " + idValue;
    }


}
