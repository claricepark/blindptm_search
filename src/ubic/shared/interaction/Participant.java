/*
 * @(#)Participant.java
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
 * @file Participant.java
 * Class file for Participant.class
 * @author John Ling
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.io.*;
import java.util.*;
import org.apache.log4j.*;


public class Participant {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Participant.class);

    // Private members of Participant class
    private String role = null;
    private boolean isTagged = false;
    private boolean isOverexpressed = false;

    private String confidenceValue = null;
    private String confidenceUnit = null;

    private int interactionId = UbicConstants.INVALIDINTID;
    
    Interactor interactor = null;

    private List features = new ArrayList();


    /**
     * @fn Participant ()
     * Constructor for Participant
     * @param interactor - An Interactor
     * @param interacionId - an interactionId
     */
    public Participant (Interactor interactor, int interactionId) {
        this.interactor = interactor;
        this.interactionId = interactionId;
    }


    /**
     * @fn Participant ()
     * Constructor for Participant
     * @param interactor - An Interactor
     */
    public Participant (Interactor interactor) {
        this.interactor = interactor;
    }


    /**
     * @fn boolean equals(Object obj)
     * Return true if this participant has the same interaction_id and 
     * interactor_id as the obj, otherwise return false.
     * @return true if this participants equals to the obj.
     */
    public boolean equals(Participant obj) {
        if (! (obj instanceof Participant)) {
            return false;
        }

        Participant other = (Participant) obj;
        return interactionId == other.getInteractionId() &&
                getInteractorId() == other.getInteractorId();

    }


    // Getters


    /**
     * @fn int getInteractionId()
     * Return the interaction_id of the Interaction that the participant involves.
     * @return The interaction_id.
     */
    public int getInteractionId() {
        return interactionId;
    }

    
    /**
     * @fn Iterator getFeatures()
     * Return an iterator of Features for this feature.
     * @return An iterator of Features.
     */
    public Iterator getFeatures() {
        return features.iterator();
    }

    /**
     * @fn String getConfidenceValue()
     * Gets the confidence_value for this participant.
     */
    public String getConfidenceValue() {
        return confidenceValue;
    }



    /**
     * @fn String getConfidenceUnit()
     * Gets the confidence_unit for this participant.
     */
    public String getConfidenceUnit() {
        return confidenceUnit;
    }


    /**
     * @fn int getInteractorId()
     * Return the interactor_id of this participant.
     * @return The interactor_id.
     */
    public int getInteractorId() {
        return interactor.getInteractorId();
    }


    /**
     * @fn String getRole ()
     * Return the role the participant plays in the interaction.
     * @return The role of a participant.
     */
    public String getRole() {
        return this.role;
    }

    /**
     * @fn boolean isTagged ()
     * Test if the participant was tagged in the experiment.
     * @return true if the participant was tagged, false otherwise
     */
    public boolean isTagged() {
        return this.isTagged;
    }

    /**
     * @fn boolean isOverexpressed ()
     * Test if the participant was overexpressed in the experiment.
     * @return true if the participant was overexpressed, false otherwise
     */
    public boolean isOverexpressed() {
        return this.isOverexpressed;
    }

    /**
     * @fn Interactor getInteractor ()
     * Return the interactor of this participant.
     * @return An interactor.
     */
    public Interactor getInteractor() {
        return this.interactor;
    }

    // Setters

    /**
     * @fn void setRole (String role)
     * Set the role of this participant.
     * @param role A role string.
     * @return void
     */
    public void setRole(String role) {
        this.role = role;
    }

    /**
     * @fn void setIsTagged (boolean isTagged)
     * Set the boolean state of the isTagged member variable.
     * @param isTagged A boolean state indicating whether the participant is tagged or not.
     * @return void
     */
    public void setIsTagged(boolean isTagged) {
        this.isTagged = isTagged;
    }

    /**
     * @fn void setIsOverexpressed (boolean isOverexpressed)
     * Set the boolean state of the isOverexpressed member variable.
     * @param isOverexpressed A boolean state of overexpression.
     * @return void
     */
    public void setIsOverexpressed(boolean isOverexpressed) {
        this.isOverexpressed = isOverexpressed;
    }

    
    /**
     * @fn void setConfidenceValue(String confidenceValue)
     * Sets the confidence value for this interaction.
     * @param confidenceValue A string for confidence value.
     */
    public void setConfidenceValue(String confidenceValue) {
        this.confidenceValue = confidenceValue;
    }



    /**
     * @fn void setConfidenceUnit(String confidenceUnit)
     * Sets the confidence_unit for this interaction.
     * @param confidenceUnit A string for confidence_unit.
     */
    public void setConfidenceUnit(String confidenceUnit) {
        this.confidenceUnit = confidenceUnit;
    }


    /**
     * @fn void setInteractor (Interactor interactor)
     * Set the interactor that is associated with this participant.
     * @param interactor An interactor.
     * @return void
     */
    public void setInteractor(Interactor interactor) {
        this.interactor = interactor;
    }
    
    
    /**
     * @fn void addFeature(Feature f)
     * Adds a Feature to the features list.
     * @param f - A Feature
     */
    public void addFeature(Feature f) {
        if (f != null) {
            features.add(f);
        }
    }

}
