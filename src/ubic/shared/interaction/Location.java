/*
 * @(#)Location.java
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
 * @file Location.java
 * Class file for Location class
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.io.*;
import java.util.*;
import org.apache.log4j.*;


public class Location {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Location.class);

    // Private members of Location class
    private int locationId = UbicConstants.INVALIDINTID;
    private int featureId = UbicConstants.INVALIDINTID;
    private int begin = UbicConstants.INVALIDINTVALUE;
    private int beginInterval = UbicConstants.INVALIDINTVALUE;
    private int end = UbicConstants.INVALIDINTVALUE;
    private int endInterval = UbicConstants.INVALIDINTVALUE;
    private int position = UbicConstants.INVALIDINTVALUE;
    private int positionInterval = UbicConstants.INVALIDINTVALUE;
    private int site = UbicConstants.INVALIDINTVALUE;


    /**
     * @fn Location ()
     * Constructor for Location
     */
    public Location () {
    }




    // Getters

    /**
     * @fn int getLocationId()
     * Return the location id for this location.
     * @return A location id.
     */
    public int getLocationId() {
        return locationId;
    }



    /**
     * @fn int getFeatureId()
     * Return the feature_id for this location.
     * @return A feature id.
     */
    public int getFeatureId() {
        return featureId;
    }



    /**
     * @fn int getBegin()
     * Return the begin for this location.
     * @return A begin.
     */
    public int getBegin() {
        return begin;
    }



    /**
     * @fn int getBeginInterval()
     * Return the beginInterval for this location.
     * @return A beginInterval.
     */
    public int getBeginInterval() {
        return beginInterval;
    }


    /**
     * @fn int getEnd()
     * Return the end for this location.
     * @return An end.
     */
    public int getEnd() {
        return end;
    }



    /**
     * @fn int getEndInterval()
     * Return the endinterval for this location.
     * @return An endinterval.
     */
    public int getEndInterval() {
        return endInterval;
    }


    /**
     * @fn int getPosition()
     * Return the position for this location.
     * @return A position.
     */
    public int getPosition() {
        return position;
    }



    /**
     * @fn int getPositionInterval()
     * Return the positioninterval for this location.
     * @return A positioninterval.
     */
    public int getPositionInterval() {
        return positionInterval;
    }

    
    /**
     * @fn int getSite()
     * Return the site for this location.
     * @return A site.
     */
    public int getSite() {
        return site;
    }



    // Setters

    /**
     * @fn void setLocationId (int locationId)
     * Sets the location id for this location.
     * @param locationId An location id.
     * @return void
     */
    public void setLocationId(int locationId) {
        this.locationId = locationId;
    }


    /**
     * @fn void setFeatureId(int featureId)
     * Sets the feature_id for this location.
     * @param featureId An feature_id.
     * @return void
     */
    public void setFeatureId(int featureId) {
        this.featureId = featureId;
    }

    /**
     * @fn void setBegin(int value)
     * Sets the begin value for this location.
     * @param value Value for begin.
     * @return void
     */
    public void setBegin(int value) {
        begin = value;
    }


    /**
     * @fn void setBeginInterval(int value)
     * Sets the beginInterval value for this location.
     * @param value Value for beginInterval.
     * @return void
     */
    public void setBeginInterval(int value) {
        beginInterval = value;
    }


    /**
     * @fn void setEndInterval(int value)
     * Sets the endinterval value for this location.
     * @param value Value for endinterval.
     * @return void
     */
    public void setEndInterval(int value) {
        endInterval = value;
    }


    /**
     * @fn void setEnd(int value)
     * Sets the end value for this location.
     * @param value Value for end.
     * @return void
     */
    public void setEnd(int value) {
        end = value;
    }


    /**
     * @fn void setPosition(int value)
     * Sets the position value for this location.
     * @param value Value for position.
     * @return void
     */
    public void setPosition(int value) {
        position = value;
    }


    /**
     * @fn void setPositionInterval(int value)
     * Sets the positioninterval value for this location.
     * @param value Value for positioninterval.
     * @return void
     */
    public void setPositionInterval(int value) {
        positionInterval = value;
    }
    
    
    /**
     * @fn void setSite(int value)
     * Sets the site value for this location.
     * @param value Value for site.
     * @return void
     */
    public void setSite(int value) {
        site = value;
    }


}
