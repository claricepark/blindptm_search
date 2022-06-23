package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class SampleHasOrganism implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String sampleId;

    /** nullable persistent field */
    private String organismId;

    /** full constructor */
    public SampleHasOrganism(String sampleId, String organismId) {
        this.sampleId = sampleId;
        this.organismId = organismId;
    }

    /** default constructor */
    public SampleHasOrganism() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSampleId() {
        return this.sampleId;
    }

    public void setSampleId(String sampleId) {
        this.sampleId = sampleId;
    }

    public String getOrganismId() {
        return this.organismId;
    }

    public void setOrganismId(String organismId) {
        this.organismId = organismId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof SampleHasOrganism) ) return false;
        SampleHasOrganism castOther = (SampleHasOrganism) other;
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
