package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Sample implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String sampleName;

    /** nullable persistent field */
    private String sampleType;

    /** nullable persistent field */
    private String donorId;

    /** nullable persistent field */
    private String donorGender;

    /** nullable persistent field */
    private String collaborator;

    /** nullable persistent field */
    private String description;

    /** full constructor */
    public Sample(String sampleName, String sampleType, String donorId, String donorGender, String collaborator, String description) {
        this.sampleName = sampleName;
        this.sampleType = sampleType;
        this.donorId = donorId;
        this.donorGender = donorGender;
        this.collaborator = collaborator;
        this.description = description;
    }

    /** default constructor */
    public Sample() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSampleName() {
        return this.sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public String getSampleType() {
        return this.sampleType;
    }

    public void setSampleType(String sampleType) {
        this.sampleType = sampleType;
    }

    public String getDonorId() {
        return this.donorId;
    }

    public void setDonorId(String donorId) {
        this.donorId = donorId;
    }

    public String getDonorGender() {
        return this.donorGender;
    }

    public void setDonorGender(String donorGender) {
        this.donorGender = donorGender;
    }

    public String getCollaborator() {
        return this.collaborator;
    }

    public void setCollaborator(String collaborator) {
        this.collaborator = collaborator;
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
        if ( !(other instanceof Sample) ) return false;
        Sample castOther = (Sample) other;
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
