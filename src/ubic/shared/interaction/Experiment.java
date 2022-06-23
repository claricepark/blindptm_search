/*
 * @(#)Experiment.java
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
 * @file Experiment.java
 * This is the source file for ubic.shared.interaction.Experiment
 *
 * @author Tao Xu, John Ling
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import java.util.*;
import ubic.shared.constants.UbicConstants;
import org.apache.log4j.*;


/**
 * @class Experiment
 * This class is models biological experiments
 *
 */
public class Experiment {

    private int experimentId = UbicConstants.INVALIDINTID;
    private String shortLabel = null;
    private String fullName = null;
    //private Organism hostOrganism = null;
    private String description = null;
    
    private String confidenceValue = null;
    private String confidenceUnit = null;
    private List dbxrefs = new ArrayList();
    private List detections = new ArrayList();
    private List attributes = new ArrayList();

    private static Logger logger = Logger.getLogger(Experiment.class);

    private Set hostOrganisms = new HashSet();

    // Setters

    /**
     * @fn void setExperimentId (int experimentId)
     * Sets the experiment id for this experiment.
     * @param experimentId A experiment identifier.
     */
    public void setExperimentId(int experimentId) {
        this.experimentId = experimentId;
    }


    /**
     * @fn void setConfidenceValue(String confidenceValue)
     * Sets the confidence_value for this experiment.
     * @param confidenceValue A string for confidence value.
     */
    public void setConfidenceValue(String confidenceValue) {
        this.confidenceValue = confidenceValue;
    }


    /**
     * @fn void setConfidenceUnit(String confidenceUnit)
     * Sets the confidence_unit for this experiment.
     * @param confidenceUnit - A string for confidence_unit.
     */
    public void setConfidenceUnit(String confidenceUnit) {
        this.confidenceUnit = confidenceUnit;
    }

    /**
     * @fn void setShortLabel (String shortLabel)
     * Sets the short label for this experiment.
     * @param shortLabel A short label string.
     */
    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }

    /**
     * @fn void setFullName (String fullName)
     * Sets the full name for this experiment.
     * @param fullName A fullname string.
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }

    /**
     * @fn void addHostOrganism(Organism host)
     * Add a host organism to this experiment.
     * @param host A host organism.
     */
    public void addHostOrganism(Organism host) {
        hostOrganisms.add(host);
    }

    /**
     * @fn void setDescription (String description)
     * Sets the description for this experiment.
     * @param description A description string.
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * @fn void addDbxref (Dbxref dbxref)
     * Adds a dbxref to the dbxrefs list.
     * @param dbxref A dbxref to be added to the dbxrefs list.
     */
    public void addDbxref(Dbxref dbxref) {
        if (dbxref != null) {
            this.dbxrefs.add(dbxref);
        }
    }

    /**
     * @fn void addDetection (Detection detection)
     * Adds a detection object to the detections list.
     * @param detection A detection object to be added to the list
     * of detections.
     */
    public void addDetection (Detection detection) {
        if (detection != null) {
            this.detections.add(detection);
        }
    }

    /**
     * @fn void addAttribute (Attribute attribute)
     * Adds an attribute object to the attributes list.
     * @param attribute A attribute object to be added to the list
     * of attributes.
     */
    public void addAttribute (Attribute attribute) {
        if (attribute != null) {
            this.attributes.add(attribute);
        }
    }

    // Getters

    /**
     * @fn int getExperimentId ()
     * Returns the experiment id for this experiment.
     */
    public int getExperimentId() {
        return this.experimentId;
    }
    
    

    /**
     * @fn String getConfidenceValue()
     * Return the confidence_value for this experiment.
     * @return The experiment's confidence_value.
     */
    public String getConfidenceValue() {
        return confidenceValue;
    }

    
    /**
     * @fn String getConfidenceUnit()
     * Return the confidence_unit for this experiment.
     * @return The experiment's confidence_unit.
     */
    public String getConfidenceUnit() {
        return confidenceUnit;
    }
    /**
     * @fn String getShortLabel ()
     * Return the short label for this experiment.
     * @return The experiment's short label.
     */
    public String getShortLabel() {
        return shortLabel;
    }

    /**
     * @fn String getFullName ()
     * Return the full name for this experiment.
     * @return The experiment's full name.
     */
    public String getFullName() {
        return fullName;
    }

    /**
     * @fn Iterator getHostOrganism ()
     * Return the host organism for this experiment.
     * @return The host organism used in this experiment.
     */
    public Iterator getHostOrganisms() {
        return hostOrganisms.iterator();
    }

    /**
     * @fn String getDescription ()
     * Return the description for this experiment.
     * @return The experiment's description.
     */
    public String getDescription() {
        return description;
    }

    /**
     * @fn Iterator getDbxrefs ()
     * Return an iterator of dbxref objects associated with this experiment.
     * @return Iterator of dbxref objects for this experiment.
     */
    public Iterator getDbxrefs() {
        return this.dbxrefs.iterator();
    }

    /**
     * @fn Iterator getDetections ()
     * Return an iterator for the detections associated with this experiment.
     * @return Iterator of detection objects for this experiment.
     */
    public Iterator getDetections() {
        return this.detections.iterator();
    }

    /**
     * @fn Iterator getAttributes ()
     * Return an iterator for the attributes associated with this experiment.
     * @return Iterator of attribute objects for this experiment.
     */
    public Iterator getAttributes () {
        return this.attributes.iterator();
    }

    /**
     * @fn void clearDbxrefs ()
     * Removes all of the dbxrefs from this experiment.
     * @exception UnsupportedOperationException
     */
    public void clearDbxrefs()
        throws UnsupportedOperationException {
        this.dbxrefs.clear();
    }

    /**
     * @fn void clearDetections ()
     * Removes all of the detections from this experiment.
     * @exception UnsupportedOperationException
     */
    public void clearDetections()
        throws UnsupportedOperationException {
        this.detections.clear();
    }

    /**
     * @fn void clearAttributes ()
     * Removes all of the attributes from this experiment.
     * @exception UnsupportedOperationException
     */
    public void clearAttributes()
        throws UnsupportedOperationException {
        this.attributes.clear();
    }
}
