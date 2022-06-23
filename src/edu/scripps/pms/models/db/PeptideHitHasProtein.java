package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class PeptideHitHasProtein implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String proteinId;

    /** nullable persistent field */
    private String peptideHitId;

    /** full constructor */
    public PeptideHitHasProtein(String proteinId, String peptideHitId) {
        this.proteinId = proteinId;
        this.peptideHitId = peptideHitId;
    }

    /** default constructor */
    public PeptideHitHasProtein() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getProteinId() {
        return this.proteinId;
    }

    public void setProteinId(String proteinId) {
        this.proteinId = proteinId;
    }

    public String getPeptideHitId() {
        return this.peptideHitId;
    }

    public void setPeptideHitId(String peptideHitId) {
        this.peptideHitId = peptideHitId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof PeptideHitHasProtein) ) return false;
        PeptideHitHasProtein castOther = (PeptideHitHasProtein) other;
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
