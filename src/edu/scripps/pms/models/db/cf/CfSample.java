package edu.scripps.pms.models.db.cf;

import java.io.Serializable;
import java.util.Date;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class CfSample implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String sampleName;

    /** nullable persistent field */
    private String treatmentId;

    /** nullable persistent field */
    private String ipProtocolId;

    /** nullable persistent field */
    private String numDishes;

    /** nullable persistent field */
    private String cellLine;

    /** nullable persistent field */
    private String preparedBy;

    /** nullable persistent field */
    private String handledBy;

    /** nullable persistent field */
    private String inputBy;

    /** nullable persistent field */
    private Date date;

    /** nullable persistent field */
    private Date modifiedTime;

    /** full constructor */
    public CfSample(String sampleName, String treatmentId, String ipProtocolId, String numDishes, String cellLine, String preparedBy, String handledBy, String inputBy, Date date, Date modifiedTime) {
        this.sampleName = sampleName;
        this.treatmentId = treatmentId;
        this.ipProtocolId = ipProtocolId;
        this.numDishes = numDishes;
        this.cellLine = cellLine;
        this.preparedBy = preparedBy;
        this.handledBy = handledBy;
        this.inputBy = inputBy;
        this.date = date;
        this.modifiedTime = modifiedTime;
    }

    /** default constructor */
    public CfSample() {
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

    public String getTreatmentId() {
        return this.treatmentId;
    }

    public void setTreatmentId(String treatmentId) {
        this.treatmentId = treatmentId;
    }

    public String getIpProtocolId() {
        return this.ipProtocolId;
    }

    public void setIpProtocolId(String ipProtocolId) {
        this.ipProtocolId = ipProtocolId;
    }

    public String getNumDishes() {
        return this.numDishes;
    }

    public void setNumDishes(String numDishes) {
        this.numDishes = numDishes;
    }

    public String getCellLine() {
        return this.cellLine;
    }

    public void setCellLine(String cellLine) {
        this.cellLine = cellLine;
    }

    public String getPreparedBy() {
        return this.preparedBy;
    }

    public void setPreparedBy(String preparedBy) {
        this.preparedBy = preparedBy;
    }

    public String getHandledBy() {
        return this.handledBy;
    }

    public void setHandledBy(String handledBy) {
        this.handledBy = handledBy;
    }

    public String getInputBy() {
        return this.inputBy;
    }

    public void setInputBy(String inputBy) {
        this.inputBy = inputBy;
    }

    public Date getDate() {
        return this.date;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public Date getModifiedTime() {
        return this.modifiedTime;
    }

    public void setModifiedTime(Date modifiedTime) {
        this.modifiedTime = modifiedTime;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof CfSample) ) return false;
        CfSample castOther = (CfSample) other;
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
