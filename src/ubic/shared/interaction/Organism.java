/*
 * @(#)Organism.java
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
 * @file Organism.java
 * This is the source file for ubic.shared.interaction.Organism
 *
 * @author Tao Xu, John Ling
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;


import ubic.shared.constants.UbicConstants;

import org.apache.log4j.*;

import java.util.*;




/**
 * @class Organism
 * This class is models an organism entity
 *
 */
public class Organism {

    private int organismId = UbicConstants.INVALIDINTID;
    private String shortLabel = null;
    private String fullName = null;
    private int taxonId = UbicConstants.INVALIDINTID;

    private Set cellTypes = new HashSet();
    private Set tissueTypes = new HashSet();
    private Set compartments = new HashSet();


    private static Logger logger = Logger.getLogger(Organism.class);


    /**
     * @fn Organism()
     * Constructor for Interactor
     */
    public Organism() {

    }

    /**
     * @fn Organism(int organismId, int taxonId, String fullName,
                                                String shortLabel)
     * Constructor for Interactor
     */
    public Organism(int organismId, int taxonId, String fullName,
                                                String shortLabel) {
        this.organismId = organismId;
        this.taxonId = taxonId;
        this.fullName = fullName;
        this.shortLabel = shortLabel;
    }

    /**
     * @fn Organism(int organismId, int taxonId, String fullName,
                                                String shortLabel)
     * Constructor for Organism
     */

    public Organism(int taxonId, String fullName, String shortLabel) {
        this.taxonId = taxonId;
        this.fullName = fullName;
        this.shortLabel = shortLabel;
    }



    /**
     * @fn void setOrganismId(int organismId)
     * Set the organism id of this Organism.
     * @return void.
     */
    public void setOrganismId(int organismId) {
        this.organismId = organismId;
    }


    /**
     * @fn void setTaxonId(int taxonId)
     * Set the taxonId of this Organism.
     * @return void.
     */
    public void setTaxonId(int taxonId) {
        this.taxonId = taxonId;
    }


    /**
     * @fn void setFullName(String fullName)
     * Set the fullName of this Organism.
     * @return void.
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }


    /**
     * @fn void setFullName(String fullName)
     * Set the fullName of this Organism.
     * @return void.
     */
    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }


    /**
     * @fn int getOrganismId()
     * Return the organism id of this Organism.
     * @return the organism id.
     */
    public int getOrganismId() {
        return organismId;
    }


    /**
     * @fn String getShortLabel()
     * Return the short label of this Organism.
     * @return the short label.
     */
    public String getShortLabel() {
        return shortLabel;
    }


    /**
     * @fn String getFullName()
     * Return the full name of this Organism.
     * @return the full name.
     */
    public String getFullName() {
        return fullName;
    }


    /**
     * @fn int getTaxonId()
     * Return the taxonomy id of this Organism.
     * @return the taxonomy id.
     */
    public int getTaxonId() {
        return taxonId;
    }


    /**
     * @fn void addCellType(String cellType)
     * Adds a cell type to the cellTypes list.
     * @param cellType - A String for a cel type.
     * @return void
     */
    public void addCellType(Localization cellType) {
        if (cellType != null) {
            cellTypes.add(cellType);
        }
    }



    /**
     * @fn void addTissueType(String tissueType)
     * Adds a tissue type to the tissueTypes list.
     * @param tissueType - A String for a tissue type.
     * @return void
     */
    public void addTissueType(Localization tissueType) {
        if (tissueType != null) {
            tissueTypes.add(tissueType);
        }
    }


    /**
     * @fn void addCompartment(String compartment)
     * Adds a compartment to the compartments list.
     * @param compartment - A String for a compartment type.
     * @return void
     */
    public void addCompartment(Localization compartment) {
        if (compartment != null) {
            compartments.add(compartment);
        }
    }

    
    /**
     * @fn Iterator getTissueTypes()
     * Return an Iterator for the tissueTypes of this Interactor.
     * @return An Iterator for the tissueTypes .
     */
    public Iterator getTissueTypes() {
        return tissueTypes.iterator();
    }

    /**
     * @fn Iterator getCellTypes()
     * Return an Iterator for the cellTypes of this Interactor.
     * @return An Iterator for the cellTypes.
     */
    public Iterator getCellTypes() {
        return cellTypes.iterator();
    }

    /**
     * @fn Iterator getCompartments()
     * Return the compartments of this Interactor.
     * @return the compartments .
     */
    public Iterator getCompartments() {
        return compartments.iterator();
    }


    /**
     * @fn boolean equals(Organism o)
     * Return true if this equals to o, otherwise return false.
     * @return true if this equals to o, otherwise false.
     * @note two Organisms are considered equal if they have the same taxonIds
     * that are not 0 or UbicConstants.INVALIDINTID. If the taxonIds are 0 or
     * UbicConstants.INVALIDINTID, then check if they have the same shortlabels
     * or fullnames.
     */
    public boolean equals(Organism o) {

        logger.debug("in equals");
        if (o == null) {
            return false;
        }

        if (taxonId != 0 && taxonId != UbicConstants.INVALIDINTID) {
            return taxonId == o.getTaxonId();
        }

        if (shortLabel != null) {
            return shortLabel.equals(o.getShortLabel());
        }

        if (fullName != null) {
            return fullName.equals(o.getFullName());
        }

        return false;
    }

    public int hashCode() {
    
        logger.debug("in hashCode");

        int hashCode = UbicConstants.INVALIDINTID;
        if (taxonId != 0 && taxonId != UbicConstants.INVALIDINTID) {
            hashCode = taxonId;
        } else {

            (shortLabel + fullName).hashCode();
        }

        return hashCode;

    }
}
