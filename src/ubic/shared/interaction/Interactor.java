/*
 * @(#)Interactor.java
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
 * @file Interactor.java
 * This is the source file for ubic.shared.interaction.Interactor
 *
 * @author Tao Xu, John Ling
 * @date $Date: 2005/05/24 22:49:29 $
 */


package ubic.shared.interaction;


import ubic.shared.constants.UbicConstants;

import org.apache.log4j.*;

import java.util.*;




/**
 * @class Interactor
 * This class is models a molecular interactor in a biological molecular
 * interaction
 *
 */
public class Interactor {

    private int interactorId = UbicConstants.INVALIDINTID;
    private String shortLabel = null;
    private String fullName = null;
    //private String organismShortLabel = null;
    //private String organismFullName = null;
    //private int taxonId = UbicConstants.INVALIDINTID;
    private String sequence = null;
    private int moltype = UbicConstants.INVALIDINTVALUE;
    private List dbxrefs = new ArrayList();
    //private Set cellTypes = new HashSet();
    //private Set tissueTypes = new HashSet();
    //private Set compartments = new HashSet();


    private Organism organism = null;

    private static Logger logger = Logger.getLogger(Interactor.class);



    /**
     * @fn Interactor (Organism org)
     * Constructor for Interactor
     * @param org Organism associated with this Interactor
     * @note org cannot be null
     */
    public Interactor (Organism org) {
        this.organism = org;
    }

    /**
     * @fn void setOrganism(Organism organism)
     * Sets the organism for this interactor.
     * @param organism - An organism.
     */
    public void setOrganism(Organism organism) {

        this.organism = organism;
    }


    /**
     * @fn void setShortLabel(String shortLabel)
     * Sets the short label for this interactor.
     * @param shortLabel - A string for short label.
     * @return void
     */
    public void setShortLabel(String shortLabel) {
        this.shortLabel = shortLabel;
    }

    /**
     * @fn void setFullName(String fullName)
     * Sets the full name for this interactor.
     * @param fullName - A string for full name.
     * @return void
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }


    /**
     * @fn void setSequence(String sequence)
     * Sets the sequence for this interactor.
     * @param sequence - A String for sequence.
     * @return void
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }


    /**
     * @fn void setMoltype(int moltype)
     * Sets the moltype for this interactor.
     * @param moltype - A int for moltype.
     * @return void
     */
    public void setMoltype(int moltype) {
        this.moltype = moltype;
    }


    /**
     * @fn void setInteractorId(int interactorId)
     * Sets the interactor id for this interactor.
     * @param interactorId - A int for interactor id.
     * @return void
     */
    public void setInteractorId(int interactorId) {
        this.interactorId = interactorId;
    }


    /**
     * @fn int getInteractorId()
     * Return the interactor id of this Interactor.
     * @return the interactor id.
     */
    public int getInteractorId() {
        return interactorId;
    }


    /**
     * @fn String getShortLabel()
     * Return the short label of this Interactor.
     * @return the short label.
     */
    public String getShortLabel() {
        return shortLabel;
    }


    /**
     * @fn String getFullName()
     * Return the full name of this Interactor.
     * @return the full name.
     */
    public String getFullName() {
        return fullName;
    }


    /**
     * @fn int getTaxonId()
     * Return the taxonomy id of this Interactor.
     * @return the taxonomy id.
     */
    public int getTaxonId() {
        return organism.getTaxonId();
    }


    /**
     * @fn int getMoltype()
     * Return the moltype of this Interactor.
     * @return the moltype.
     */
    public int getMoltype() {
        return moltype;
    }


    /**
     * @fn String getSequence()
     * Return the sequence of this Interactor.
     * @return the sequence.
     */
    public String getSequence() {
        return sequence;
    }


    /**
     * @fn String organismShortLabel()
     * Return the organism short label of this Interactor.
     * @return the organism short label.
     */
    public String getOrganismShortLabel() {
        return organism.getShortLabel();
    }


    /**
     * @fn String organismFullName()
     * Return the organism full name of this Interactor.
     * @return the organism full name.
     */
    public String getOrganismFullName() {
        return organism.getFullName();
    }

    /**
     * @fn Iterator getDbxrefs()
     * Return An Iterator for the dbxrefs of this Interactor.
     * @return An Iterator for the dbxrefs.
     */
    public Iterator getDbxrefs() {
        return dbxrefs.iterator();
    }


    /**
     * @fn Iterator getTissueTypes()
     * Return an Iterator of Localizations for tissueTypes of this Interactor.
     * @return An Iterator for the tissueTypes .
     */
    public Iterator getTissueTypes() {
        return organism.getTissueTypes();
    }

    /**
     * @fn Iterator getCellTypes()
     * Return an Iterator of Localizations for cellTypes of this Interactor.
     * @return An Iterator for the cellTypes.
     */
    public Iterator getCellTypes() {
        return organism.getCellTypes();
    }

    /**
     * @fn Iterator getCompartments()
     * Return an Iterator of Localizations for compartments of this Interactor.
     * @return the compartments .
     */
    public Iterator getCompartments() {
        return organism.getCompartments();
    }

    /**
     * @fn Iterator getDbxref(String idType, String idSource)
     * Return the Dbxref object that matches the given idType and idSource.
     * This makes the assumption that there is only one unique dbxref in
     * this Interactor's list of dbxref's, unique in the sense of idType and
     * idSource.
     * @return A dbxref object.
     */
    public Dbxref getDbxref(String idType, String dbSource) {
        Iterator dbxrefIterator = this.getDbxrefs();

        Dbxref d = null;
        while (dbxrefIterator.hasNext()) {
            d = (Dbxref)dbxrefIterator.next();

            if (d.getIdType().toUpperCase().equals(idType.toUpperCase()) && d.getDbSource().toUpperCase().equals(dbSource.toUpperCase())) {
                return d;
            }
        }

        // Return null when no such dbxref was found.
        return null;
    }


    /**
     * @fn Organism getOrganism()
     * Return the organism of this interactor.
     * @return An Organism
     */
    public Organism getOrganism() {
        return organism;
    }

    /**
     * @fn void addDbxref(Dbxref dbxref)
     * Adds a dbxref to the dbxrefs list.
     * @param dbxref - A Dbxref.
     * @return void
     */
    public void addDbxref(Dbxref dbxref) {
        if (dbxref != null) {
            dbxrefs.add(dbxref);
        }
    }


    /**
     * @fn void addCellType(String cellType)
     * Adds a cell type to the cellTypes list.
     * @param cellType - A Localization object for a cell type.
     * @return void
     */
    public void addCellType(Localization cellType) {
        if (cellType != null) {
            organism.addCellType(cellType);
        }
    }



    /**
     * @fn void addTissueType(String tissueType)
     * Adds a tissue type to the tissueTypes list.
     * @param tissueType - A Localization object for a tissue type.
     * @return void
     */
    public void addTissueType(Localization tissueType) {
        if (tissueType != null) {
            organism.addTissueType(tissueType);
        }
    }


    /**
     * @fn void addCompartment(String compartment)
     * Adds a compartment to the compartments list.
     * @param compartment - A Localization object for a compartment type.
     * @return void
     */
    public void addCompartment(Localization compartment) {
        if (compartment != null) {
            organism.addCompartment(compartment);
        }
    }
}
