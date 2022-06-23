package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Ptm implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private Integer position;

    /** nullable persistent field */
    private String peptideHitId;

    /** nullable persistent field */
    private String modificationTypeId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.ModificationType modificationType;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.PeptideHit peptideHit;

    /** full constructor */
    public Ptm(Integer position, String peptideHitId, String modificationTypeId, edu.scripps.pms.models.db.ModificationType modificationType, edu.scripps.pms.models.db.PeptideHit peptideHit) {
        this.position = position;
        this.peptideHitId = peptideHitId;
        this.modificationTypeId = modificationTypeId;
        this.modificationType = modificationType;
        this.peptideHit = peptideHit;
    }

    /** default constructor */
    public Ptm() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public Integer getPosition() {
        return this.position;
    }

    public void setPosition(Integer position) {
        this.position = position;
    }

    public String getPeptideHitId() {
        return this.peptideHitId;
    }

    public void setPeptideHitId(String peptideHitId) {
        this.peptideHitId = peptideHitId;
    }

    public String getModificationTypeId() {
        return this.modificationTypeId;
    }

    public void setModificationTypeId(String modificationTypeId) {
        this.modificationTypeId = modificationTypeId;
    }

    public edu.scripps.pms.models.db.ModificationType getModificationType() {
        return this.modificationType;
    }

    public void setModificationType(edu.scripps.pms.models.db.ModificationType modificationType) {
        this.modificationType = modificationType;
    }

    public edu.scripps.pms.models.db.PeptideHit getPeptideHit() {
        return this.peptideHit;
    }

    public void setPeptideHit(edu.scripps.pms.models.db.PeptideHit peptideHit) {
        this.peptideHit = peptideHit;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Ptm) ) return false;
        Ptm castOther = (Ptm) other;
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
