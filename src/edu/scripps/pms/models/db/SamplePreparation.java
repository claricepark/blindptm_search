package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class SamplePreparation implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String protocol;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.SeparationType separationTypeId;

    /** full constructor */
    public SamplePreparation(String protocol, edu.scripps.pms.models.db.SeparationType separationTypeId) {
        this.protocol = protocol;
        this.separationTypeId = separationTypeId;
    }

    /** default constructor */
    public SamplePreparation() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getProtocol() {
        return this.protocol;
    }

    public void setProtocol(String protocol) {
        this.protocol = protocol;
    }

    public edu.scripps.pms.models.db.SeparationType getSeparationTypeId() {
        return this.separationTypeId;
    }

    public void setSeparationTypeId(edu.scripps.pms.models.db.SeparationType separationTypeId) {
        this.separationTypeId = separationTypeId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof SamplePreparation) ) return false;
        SamplePreparation castOther = (SamplePreparation) other;
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
