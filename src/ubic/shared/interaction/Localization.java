/*
 * @(#)Localization.java
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
 * @file Localization.java
 * Class file for Localization class
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.io.*;
import java.util.*;
import org.apache.log4j.*;


public class Localization {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Localization.class);

    // Private members of Localization class
    private int localizationId = UbicConstants.INVALIDINTID;
    private String localization = null;
    //private int dbxrefId = UbicConstants.INVALIDINTID;
    //private String xrefValue = null;

    private Dbxref dbxref = null;


    /**
     * @fn Localization ()
     * Constructor for Localization
     */
    public Localization () {
    }


    /**
     * @fn Localization ()
     * Constructor for Localization
     */
    public Localization (String localization) {
        this.localization = localization;
    }


    /**
     * @fn Localization ()
     * Constructor for Localization
     */
    public Localization (String localization, Dbxref dbxref) {
        this.localization = localization;
        this.dbxref = dbxref;
    }


    /**
     * @fn Localization ()
     * Constructor for Localization
     */
    public Localization (int locId, String localization, Dbxref dbxref) {
        localizationId = locId;
        this.localization = localization;
        this.dbxref = dbxref;
    }


    public boolean equals(Localization obj) {

        boolean isEqual = false;

        if (localization != null) {
            return localization.equals(obj.getLocalization());
        }

        return isEqual;
    }
    // Getters

    /**
     * @fn int getLocalizationId()
     * Return the localization_id for this localization.
     * @return A localization_id.
     */
    public int getLocalizationId() {
        return localizationId;
    }

    
    /**
     * @fn Dbxref getDbxref()
     * Return the Dbxref of this localization.
     * @return A Dbxref.
     */
    public Dbxref getDbxref() {
        return dbxref;
    }

    /**
     * @fn String getShortLabel ()
     * Return the localization this Localization.
     * @return The localization.
     */
    public String getLocalization() {
        return localization;
    }


    // Setters

    /**
     * @fn void setLocalizationId (int localizationId)
     * Sets the localization id for this localization.
     * @param localizationId An localization id.
     * @return void
     */
    public void setLocalizationId(int localizationId) {
        this.localizationId = localizationId;
    }


    /**
     * @fn void setLocalization(String localization)
     * Sets the localization of this Localization.
     * @param localization A string for localization.
     */
    public void setLocalization(String localization) {
        this.localization = localization;
    }


    /**
     * @fn void setDbxref(Dbxref dbxref)
     * Sets the dbxref of this Localization.
     * @param dbxref A dbxref.
     */
    public void setDbxref(Dbxref dbxref) {
        this.dbxref = dbxref;
    }


}
