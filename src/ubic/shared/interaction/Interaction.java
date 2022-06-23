/*
 * @(#)Interaction.java
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
 * @file Interaction.java
 * Class file for Interaction class
 * @author John Ling and Tao Xu
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;

import ubic.shared.constants.UbicConstants;

// java util
import java.io.*;
import java.util.*;
import org.apache.log4j.*;


public class Interaction {


    // Create a log4j logger for debugging
    private static Logger logger = Logger.getLogger(Interaction.class);

    // Private members of Interaction class
    private int interactionId = UbicConstants.INVALIDINTID;
    private String interactionType = null;
    private String shortLabel = null;
    private String fullName = null;
    private String description = null;
    private String confidenceValue = null;
    private String confidenceUnit = null;

    private List participants = new ArrayList();
    private List experiments = new ArrayList();
    private List dbxrefs = new ArrayList();


    /**
     * @fn Interactor ()
     * Constructor for Interactor
     */
    public Interaction () {
    }


    // Getters

    /**
     * @fn int getInteractor ()
     * Return the interaction id for this interaction.
     * @return An interaction id.
     */
    public int getInteractionId() {
        return this.interactionId;
    }

    /**
     * @fn String getInteractionType()
     * Return the interaction type for this interaction.
     * @return An interaction type.
     */
    public String getInteractionType() {
        return this.interactionType;
    }

    /**
     * @fn String getShortLabel ()
     * Return the short label for this interaction.
     * @return The interaction's short label.
     */
    public String getShortLabel() {
        return this.shortLabel;
    }

    /**
     * @fn String getFullName ()
     * Return the full name for this interaction.
     * @return The interaction's full name.
     */
    public String getFullName() {
        return this.fullName;
    }

    /**
     * @fn String getDescription ()
     * Return the description for this interaction.
     * @return The interaction's description.
     */
    public String getDescription() {
        return this.description;
    }

    /**
     * @fn Iterator getParticpants ()
     * Return an iterator of interacting participants involved in this interaction.
     * @return An iterator of participants.
     */
    public Iterator getParticipants() {
        return this.participants.iterator();
    }

    /**
     * @fn int getParticipantsSize ()
     * Return the number of participants in this interaction.
     * @return Number of elements in the participants list.
     */
    public int getParticipantsSize() {
        return this.participants.size();
    }

    
    
    /**
     * @fn String getConfidenceValue()
     * Gets the confidence_value for this interaction.
     */
    public String getConfidenceValue() {
        return confidenceValue;
    }



    /**
     * @fn String getConfidenceUnit()
     * Gets the confidence_unit for this interaction.
     */
    public String getConfidenceUnit() {
        return confidenceUnit;
    }

    /**
     * @fn boolean isParticipantsEmpty ()
     * Tests if the participants list in this interaction is empty or not.
     * @return true if the participants list is empty, false otherwise
     */
    public boolean isParticipantsEmpty() {
        return this.participants.isEmpty();
    }

    /**
     * @fn boolean containsParticpant (Participant participant)
     * Tests if a given participant exists in this interaction.
     * @param A participant
     * @return true if the participant exists in the interaction, false otherwise
     */
    public boolean containsParticipant(Participant participant) {
        return this.participants.contains(participant);
    }

    /**
     * @fn Iterator getExperiments ()
     * Return an iterator of experiments involved in this interaction.
     * @return An iterator of experiments.
     */
    public Iterator getExperiments() {
        return this.experiments.iterator();
    }



    /**
     * @fn Iterator getDbxrefs()
     * Return an iterator of Dbxrefs for this interaction.
     * @return An iterator of Dbxref.
     */
    public Iterator getDbxrefs() {
        return dbxrefs.iterator();
    }
    /**
     * @fn int getExperimentsSize ()
     * Return the number of experiments in this interaction.
     * @return Number of elements in the experiments list.
     */
    public int getExperimentsSize() {
        return this.experiments.size();
    }

    /**
     * @fn boolean isExperimentsEmpty ()
     * Tests if the experiments list in this interaction is empty or not.
     * @return true if the experiments list is empty, false otherwise
     */
    public boolean isExperimentsEmpty() {
        return this.experiments.isEmpty();
    }

    /**
     * @fn boolean containsExperiment (Experiment experiment)
     * Tests if a given experiment exists in this interaction.
     * @param experiment - An experiment
     * @return true if the experiment exists in the interaction, false otherwise
     */
    public boolean containsExperiment(Experiment experiment) {
        return this.experiments.contains(experiment);
    }

    // Setters

    /**
     * @fn void setInteractionId (int interactionId)
     * Sets the interaction id for this interaction.
     * @param interactionId An interaction id.
     * @return void
     */
    public void setInteractionId(int interactionId) {
        this.interactionId = interactionId;
    }

    /**
     * @fn void setInteractionType(String interactionType)
     * Sets the interaction type for this interaction.
     * @param interactionType An interaction type.
     * @return void
     */
    public void setInteractionType(String interactionType) {
        this.interactionType = interactionType;
    }

    /**
     * @fn void setShortLabel (String shortLabel)
     * Sets the short label for this interaction.
     * @param shortLabel A short label string.
     */
    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }

    /**
     * @fn void setFullName (String fullName)
     * Sets the full name for this interaction.
     * @param fullName A fullname string.
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }

    /**
     * @fn void setDescription (String description)
     * Sets the description for this interaction.
     * @param description A description string.
     */
    public void setDescription(String description) {
        this.description = description;
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
     * @fn void addParticipant (Participant participant)
     * Adds a participant to the participants list.
     * @param participant A participant.
     */
    public void addParticipant(Participant participant) {
        if (participant != null) {
            participants.add(participant);
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

    /**
     * @fn void clearParticipants ()
     * Removes all of the participants from this interaction.
     * @exception UnsupportedOperationException
     */
    public void clearParticipants()
        throws UnsupportedOperationException {
        this.participants.clear();
    }


    /**
     * @fn void addExperiment (Experiment experiment)
     * Adds an experiment to the experiments list.
     * @param experiment An experiment.
     */
    public void addExperiment(Experiment experiment) {
        if (experiment != null) {
            this.experiments.add(experiment);
        }
    }

    /**
     * @fn void clearExperiments ()
     * Removes all of the experiments from this interaction.
     * @exception UnsupportedOperationException
     */
    public void clearExperiments()
        throws UnsupportedOperationException {
        this.experiments.clear();
    }
}
