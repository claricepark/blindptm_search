/*
 * @(#)Feature.java
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
 * @file Feature.java
 * Class file for Feature class
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.io.*;
import java.util.*;
import org.apache.log4j.*;


public class Feature {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Feature.class);

    // Private members of Feature class
    private int featureId = UbicConstants.INVALIDINTID;
    private String shortLabel = null;
    private String fullName = null;

    private List dbxrefs = new ArrayList();
    private List detections = new ArrayList();
    private List locations = new ArrayList();



    /**
     * @fn Feature ()
     * Constructor for Feature
     */
    public Feature () {
    }


    /**
     * @fn Feature ()
     * Constructor for Feature
     */
    public Feature (String fullName, String shortLabel) {
        this.fullName = fullName;
        this.shortLabel = shortLabel;
    }


        /**
     * @fn Feature ()
     * Constructor for Feature
     */
    public Feature (int fid, String fullName, String shortLabel) {
        featureId = fid;
        this.fullName = fullName;
        this.shortLabel = shortLabel;
    }


    // Getters

    /**
     * @fn int getInteractor ()
     * Return the feature id for this feature.
     * @return An feature id.
     */
    public int getFeatureId() {
        return featureId;
    }

    /**
     * @fn String getShortLabel ()
     * Return the short label for this feature.
     * @return The feature's short label.
     */
    public String getShortLabel() {
        return shortLabel;
    }

    /**
     * @fn String getFullName ()
     * Return the full name for this feature.
     * @return The feature's full name.
     */
    public String getFullName() {
        return fullName;
    }

    /**
     * @fn Iterator getDbxrefs()
     * Return an iterator of Dbxrefs for this feature.
     * @return An iterator of Dbxref.
     */
    public Iterator getDbxrefs() {
        return dbxrefs.iterator();
    }


    /**
     * @fn Iterator getLocations()
     * Return an iterator of locations for this feature.
     * @return An iterator of locations.
     */
    public Iterator getLocations() {
        return locations.iterator();
    }


    /**
     * @fn Iterator getDetections()
     * Return an iterator of Detections for this feature.
     * @return An iterator of detections.
     */
    public Iterator getDetections() {
        return detections.iterator();
    }


    // Setters

    /**
     * @fn void setFeatureId (int featureId)
     * Sets the feature id for this feature.
     * @param featureId An feature id.
     * @return void
     */
    public void setFeatureId(int featureId) {
        this.featureId = featureId;
    }


    /**
     * @fn void setShortLabel (String shortLabel)
     * Sets the short label for this feature.
     * @param shortLabel A short label string.
     */
    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }

    /**
     * @fn void setFullName (String fullName)
     * Sets the full name for this feature.
     * @param fullName A fullname string.
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }



    /**
     * @fn void addLocation(Location loc)
     * Adds a Location to the locations list.
     * @param loc - A Location.
     */
    public void addLocation(Location loc) {
        if (loc != null) {
            locations.add(loc);
        }
    }

    
    /**
     * @fn void addDetection(Detection d)
     * Adds a Detection to the detections list.
     * @param d - A Detection.
     */
    public void addDetection(Detection d) {
        if (d != null) {
            detections.add(d);
        }
    }

    /**
     * @fn void addDbxref(Dbxref dbxref)
     * Adds a Dbxref to the dbxrefs list.
     * @param dbxref - A Dbxref.
     */
    public void addDbxref(Dbxref dbxref) {
        if (dbxref != null) {
            dbxrefs.add(dbxref);
        }
    }


}
