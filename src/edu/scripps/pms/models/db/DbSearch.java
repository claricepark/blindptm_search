package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Date;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class DbSearch implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String mspExperimentId;

    /** nullable persistent field */
    private String dbSearchParametersId;

    /** nullable persistent field */
    private String userId;

    /** nullable persistent field */
    private Date idDate;

    /** nullable persistent field */
    private String resultPath;

    /** full constructor */
    public DbSearch(String mspExperimentId, String dbSearchParametersId, String userId, Date idDate, String resultPath) {
        this.mspExperimentId = mspExperimentId;
        this.dbSearchParametersId = dbSearchParametersId;
        this.userId = userId;
        this.idDate = idDate;
        this.resultPath = resultPath;
    }

    /** default constructor */
    public DbSearch() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getMspExperimentId() {
        return this.mspExperimentId;
    }

    public void setMspExperimentId(String mspExperimentId) {
        this.mspExperimentId = mspExperimentId;
    }

    public String getDbSearchParametersId() {
        return this.dbSearchParametersId;
    }

    public void setDbSearchParametersId(String dbSearchParametersId) {
        this.dbSearchParametersId = dbSearchParametersId;
    }

    public String getUserId() {
        return this.userId;
    }

    public void setUserId(String userId) {
        this.userId = userId;
    }

    public Date getIdDate() {
        return this.idDate;
    }

    public void setIdDate(Date idDate) {
        this.idDate = idDate;
    }

    public String getResultPath() {
        return this.resultPath;
    }

    public void setResultPath(String resultPath) {
        this.resultPath = resultPath;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof DbSearch) ) return false;
        DbSearch castOther = (DbSearch) other;
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
