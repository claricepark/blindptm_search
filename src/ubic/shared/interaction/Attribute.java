/*
 * @(#)Attribute.java
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
 * @file Attribute.java
 * This is the source file for ubic.shared.interaction.Attribute
 *
 * @author Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;


import ubic.shared.constants.UbicConstants;

import org.apache.log4j.*;

import java.util.*;




/**
 * @class Attribute
 * This class is models a molecular attribute in a biological molecular
 * interaction
 *
 */
public class Attribute {

    private int attributeId = UbicConstants.INVALIDINTID;
    private String name = null;
    private String value = null;
    private int experimentId = UbicConstants.INVALIDINTID;



    private static Logger logger = Logger.getLogger(Attribute.class);


    

    /**
     * @fn Attribute(String name, String value)
     * Constructing an Attribute object
     * @param name - the name of this Attribute
     * @param value - the value of this Attribute
     */
    public Attribute(String name, String value) {
        this.name = name;
        this.value = value;
    }



    /**
     * @fn Attribute(int attributeId, int experimentId, String name, String value)
     * Constructing an Attribute object
     * @param attributeId - the attribute_id of this Attribute
     * @param experimentId - the experiment_id of this Attribute
     * @param name - the name of this Attribute
     * @param value - the value of this Attribute
     */
    public Attribute(int attributeId, int experimentId, String name, String value) {
        this.name = name;
        this.value = value;
        this.attributeId = attributeId;
        this.experimentId = experimentId;
    }

    /**
     * @fn void setAttributeId(int attributeId)
     * Sets the attribute id for this attribute.
     * @param attributeId - A int for attribute id.
     * @return void
     */
    public void setAttributeId(int attributeId) {
        this.attributeId = attributeId;
    }


    /**
     * @fn void setExperimentId(int experimentId)
     * Sets the experiment_id for this Attribute.
     * @param experimentId - A int for experiment_id.
     * @return void
     */
    public void setExperimentId(int experimentId) {
        this.experimentId = experimentId;
    }



    /**
     * @fn int getAttributeId()
     * Return the attribute_id of this Attribute.
     * @return the attribute id.
     */
    public int getAttributeId() {
        return attributeId;
    }


    /**
     * @fn String getName()
     * Return the name of this Attribute.
     * @return the name this this Attribute.
     */
    public String getName() {
        return name;
    }


    /**
     * @fn String getValue()
     * Return the value of this Attribute.
     * @return the value of this Attribute.
     */
    public String getValue() {
        return value;
    }


    /**
     * @fn int getExperimentId()
     * Return the experiment_id of this Attribute.
     * @return the experiment_id.
     */
    public int getExperimentId() {
        return experimentId;
    }


}
