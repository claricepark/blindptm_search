/*
 * @(#)Detection.java
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
 * @file Detection.java
 * Class file for Detection class
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.util.*;
import org.apache.log4j.*;


public class Detection {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Detection.class);

    // Private members of Detection class
    private int detectionId = UbicConstants.INVALIDINTID;
    private String shortLabel = null;
    private String fullName = null;
    private String description = null;
    private String detectionType = null;

    private List dbxrefs = new ArrayList();



    /**
     * @fn Detection (String shortLabel, String fullName, String description)
     * Constructor for Interactor
     */
    public Detection (String shortLabel, String fullName, String description) {
        this.shortLabel = shortLabel;
        this.fullName = fullName;
        this.description = description;
    }

    /**
     * @fn Detection (String shortLabel, String fullName,
        String description, String detection)
     * Constructor for Interactor
     */
    public Detection (String shortLabel, String fullName,
        String description, String detectionType) {
        this.shortLabel = shortLabel;
        this.fullName = fullName;
        this.description = description;
        this.detectionType = detectionType;
    }
    
    /**
     * @fn Detection (String shortLabel, String fullName,
        String description, String detection)
     * Constructor for Interactor
     */
    public Detection (int detectionId, String shortLabel, String fullName,
        String description, String detectionType) {
        this.shortLabel = shortLabel;
        this.fullName = fullName;
        this.description = description;
        this.detectionType = detectionType;
        this.detectionId = detectionId;
    }

    /**
     * @fn void setDetectionId(int detectionId)
     * Set the detection_id for this Detection.
     * @param detectionId - a detection_id.
     * @return void.
     */
    public void setDetectionId(int detectionId) {
        this.detectionId = detectionId;
    }


    /**
     * @fn void setDetectionType(String detectionType)
     * Set the detectionType for this Detection.
     * @param detectionType - a string for detectionType.
     * @return void.
     */
    public void setDetectionType(String detectionType) {
        this.detectionType = detectionType;
    }

    /**
     * @fn void addDbxref(Dbxref dbxref)
     * Add a Dbxref to this Detection.
     * @param dbxref - a dbxref
     * @return void.
     */
    public void addDbxref(Dbxref dbxref) {
        dbxrefs.add(dbxref);
    }

    /**
     * @fn void setShortLabel(String shortLabel)
     * Set the short label for this Detection.
     * @param shortLabel - a String for the description
     * @return void.
     */
    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }

    /**
     * @fn void setFullName(String fullName)
     * Set the full name for this Detection.
     * @param fullName - a String for the description
     * @return void.
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }

    /**
     * @fn void setDescription(String description)
     * Set the description for this Detection.
     * @param description - a String for the description
     * @return void.
     */
    public void setDescription(String description) {
        this.description = description;
    }

    // Getters

    /**
     * @fn int getDetectionId()
     * Return the detection id for this Detection.
     * @return A detection id.
     */
    public int getDetectionId() {
        return this.detectionId;
    }



    /**
     * @fn String getDetection()
     * Return the detection for this Detection.
     * @return A String for detection.
     */
    public String getDetectionType() {
        return detectionType;
    }

    /**
     * @fn List getDbxrefs()
     * Return an Iteractor for the dbxrefs of this Detection.
     * @return An Iteractor for the dbxrefs.
     */
    public Iterator getDbxrefs() {
        return dbxrefs.iterator();
    }


    /**
     * @fn String getShortLabel()
     * Return the short label of this Detection.
     * @return A String for detection short label.
     */
    public String getShortLabel() {
        return shortLabel;
    }

    /**
     * @fn String getFullName()
     * Return the full name of this Detection.
     * @return A String for detection full name.
     */
    public String getFullName() {
        return fullName;
    }

    /**
     * @fn String getDescription()
     * Return the description of this Detection.
     * @return A String for the description.
     */
    public String getDescription() {
        return description;
    }
}
