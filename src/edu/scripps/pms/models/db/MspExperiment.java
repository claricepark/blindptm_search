package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Date;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class MspExperiment implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String description;

    /** nullable persistent field */
    private Date startDate;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.Project projectId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.Users usersId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.Sample sampleId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.MassSpecInstrument massSpecInstrumentId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.SamplePreparation samplePreparationId;

    /** full constructor */
    public MspExperiment(String description, Date startDate, edu.scripps.pms.models.db.Project projectId, edu.scripps.pms.models.db.Users usersId, edu.scripps.pms.models.db.Sample sampleId, edu.scripps.pms.models.db.MassSpecInstrument massSpecInstrumentId, edu.scripps.pms.models.db.SamplePreparation samplePreparationId) {
        this.description = description;
        this.startDate = startDate;
        this.projectId = projectId;
        this.usersId = usersId;
        this.sampleId = sampleId;
        this.massSpecInstrumentId = massSpecInstrumentId;
        this.samplePreparationId = samplePreparationId;
    }

    /** default constructor */
    public MspExperiment() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Date getStartDate() {
        return this.startDate;
    }

    public void setStartDate(Date startDate) {
        this.startDate = startDate;
    }

    public edu.scripps.pms.models.db.Project getProjectId() {
        return this.projectId;
    }

    public void setProjectId(edu.scripps.pms.models.db.Project projectId) {
        this.projectId = projectId;
    }

    public edu.scripps.pms.models.db.Users getUsersId() {
        return this.usersId;
    }

    public void setUsersId(edu.scripps.pms.models.db.Users usersId) {
        this.usersId = usersId;
    }

    public edu.scripps.pms.models.db.Sample getSampleId() {
        return this.sampleId;
    }

    public void setSampleId(edu.scripps.pms.models.db.Sample sampleId) {
        this.sampleId = sampleId;
    }

    public edu.scripps.pms.models.db.MassSpecInstrument getMassSpecInstrumentId() {
        return this.massSpecInstrumentId;
    }

    public void setMassSpecInstrumentId(edu.scripps.pms.models.db.MassSpecInstrument massSpecInstrumentId) {
        this.massSpecInstrumentId = massSpecInstrumentId;
    }

    public edu.scripps.pms.models.db.SamplePreparation getSamplePreparationId() {
        return this.samplePreparationId;
    }

    public void setSamplePreparationId(edu.scripps.pms.models.db.SamplePreparation samplePreparationId) {
        this.samplePreparationId = samplePreparationId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof MspExperiment) ) return false;
        MspExperiment castOther = (MspExperiment) other;
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
