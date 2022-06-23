package edu.scripps.pms.models.db.cf;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class GeneRegulationType implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String geneRegulationType;

    /** nullable persistent field */
    private String description;

    /** full constructor */
    public GeneRegulationType(String geneRegulationType, String description) {
        this.geneRegulationType = geneRegulationType;
        this.description = description;
    }

    /** default constructor */
    public GeneRegulationType() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getGeneRegulationType() {
        return this.geneRegulationType;
    }

    public void setGeneRegulationType(String geneRegulationType) {
        this.geneRegulationType = geneRegulationType;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof GeneRegulationType) ) return false;
        GeneRegulationType castOther = (GeneRegulationType) other;
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
