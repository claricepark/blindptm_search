package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class DataSource implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String mspExperimentId;

    /** nullable persistent field */
    private String userId;

    /** nullable persistent field */
    private String rowFileName;

    /** nullable persistent field */
    private String extractedFileName;

    /** nullable persistent field */
    private String extractor;

    /** nullable persistent field */
    private String extractorVersion;

    /** full constructor */
    public DataSource(String mspExperimentId, String userId, String rowFileName, String extractedFileName, String extractor, String extractorVersion) {
        this.mspExperimentId = mspExperimentId;
        this.userId = userId;
        this.rowFileName = rowFileName;
        this.extractedFileName = extractedFileName;
        this.extractor = extractor;
        this.extractorVersion = extractorVersion;
    }

    /** default constructor */
    public DataSource() {
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

    public String getUserId() {
        return this.userId;
    }

    public void setUserId(String userId) {
        this.userId = userId;
    }

    public String getRowFileName() {
        return this.rowFileName;
    }

    public void setRowFileName(String rowFileName) {
        this.rowFileName = rowFileName;
    }

    public String getExtractedFileName() {
        return this.extractedFileName;
    }

    public void setExtractedFileName(String extractedFileName) {
        this.extractedFileName = extractedFileName;
    }

    public String getExtractor() {
        return this.extractor;
    }

    public void setExtractor(String extractor) {
        this.extractor = extractor;
    }

    public String getExtractorVersion() {
        return this.extractorVersion;
    }

    public void setExtractorVersion(String extractorVersion) {
        this.extractorVersion = extractorVersion;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof DataSource) ) return false;
        DataSource castOther = (DataSource) other;
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
