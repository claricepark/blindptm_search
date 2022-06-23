package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class ModificationType implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String modificationType;

    /** nullable persistent field */
    private String structure;

    /** nullable persistent field */
    private String massShift;

    /** nullable persistent field */
    private String description;

    /** nullable persistent field */
    private Character modifiedResidue;

    /** full constructor */
    public ModificationType(String modificationType, String structure, String massShift, String description, Character modifiedResidue) {
        this.modificationType = modificationType;
        this.structure = structure;
        this.massShift = massShift;
        this.description = description;
        this.modifiedResidue = modifiedResidue;
    }

    /** default constructor */
    public ModificationType() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getModificationType() {
        return this.modificationType;
    }

    public void setModificationType(String modificationType) {
        this.modificationType = modificationType;
    }

    public String getStructure() {
        return this.structure;
    }

    public void setStructure(String structure) {
        this.structure = structure;
    }

    public String getMassShift() {
        return this.massShift;
    }

    public void setMassShift(String massShift) {
        this.massShift = massShift;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Character getModifiedResidue() {
        return this.modifiedResidue;
    }

    public void setModifiedResidue(Character modifiedResidue) {
        this.modifiedResidue = modifiedResidue;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof ModificationType) ) return false;
        ModificationType castOther = (ModificationType) other;
        return new EqualsBuilder()
            .append(this.getId(), castOther.getId())
            .isEquals();
    }

    public int hashCode() {
        return new HashCodeBuilder()
            .append(getId())
            .toHashCode();
    }

}
