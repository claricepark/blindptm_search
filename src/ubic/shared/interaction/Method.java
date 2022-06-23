/*
 * @(#)Method.java
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
 * @file Method.java
 * Class file for Method class
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.util.*;
import org.apache.log4j.*;


public class Method {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Method.class);

    // Private members of Method class
    private int methodId = UbicConstants.INVALIDINTID;
    private String methodName = null;
    private String description = null;
    private String detection = null;

    private List dbxrefs = new ArrayList();



    /**
     * @fn Method ()
     * Constructor for Interactor
     */
    public Method (String methodName, String description) {
        this.methodName = methodName;
        this.description = description;
    }


    /**
     * @fn Method ()
     * Constructor for Interactor
     */
    public Method (String methodName, String description, String detection) {
        this.methodName = methodName;
        this.description = description;
        this.detection = detection;
    }

    
    
    /**
     * @fn Method ()
     * Constructor for Interactor
     */
    public Method (int id, String methodName, String description, String detection) {
        this.methodName = methodName;
        this.description = description;
        this.detection = detection;
        methodId = id;
    }


    /**
     * @fn void setMethodId(int methodId)
     * Set the method_id for this Method.
     * @param methodId - a method_id.
     * @return void.
     */
    public void setMethodId(int methodId) {
        this.methodId = methodId;
    }


    /**
     * @fn void setDetection(String detection)
     * Set the detection for this Method.
     * @param detection - a string fir detection.
     * @return void.
     */
    public void setDetection(String detection) {
        this.detection = detection;
    }

    /**
     * @fn void addDbxref(Dbxref dbxref)
     * Add a Dbxref to this Method.
     * @param dbxref - a dbxref
     * @return void.
     */
    public void addDbxref(Dbxref dbxref) {
        dbxrefs.add(dbxref);
    }


    /**
     * @fn void setDescription(String description)
     * Set the description for this Method.
     * @param description - a String for the description
     * @return void.
     */
    public void setDescription(String description) {
        this.description = description;
    }

    // Getters

    /**
     * @fn int getMethodId()
     * Return the method id for this Method.
     * @return A method id.
     */
    public int getMethodId() {
        return this.methodId;
    }


    
    /**
     * @fn String getDetection()
     * Return the detection for this Method.
     * @return A String for detection.
     */
    public String getDetection() {
        return detection;
    }

    /**
     * @fn List getDbxrefs()
     * Return an Iteractor for the dbxrefs of this Method.
     * @return An Iteractor for the dbxrefs.
     */
    public Iterator getDbxrefs() {
        return dbxrefs.iterator();
    }


    /**
     * @fn String getMethodName()
     * Return the name of this Method.
     * @return A String for method name.
     */
    public String getMethodName() {
        return methodName;
    }


    /**
     * @fn String getDescription()
     * Return the description of this Method.
     * @return A String for the description.
     */
    public String getDescription() {
        return description;
    }

}
