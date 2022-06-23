package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class DTASelectBestPeptide implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String peptideHitId;

    /** nullable persistent field */
    private String proteinHitId;

    /** nullable persistent field */
    private String spectraCount;

    /** nullable persistent field */
    private Boolean isUnique;

    /** nullable persistent field */
    private String confidence;

    /** full constructor */
    public DTASelectBestPeptide(String peptideHitId, String proteinHitId, String spectraCount, Boolean isUnique, String confidence) {
        this.peptideHitId = peptideHitId;
        this.proteinHitId = proteinHitId;
        this.spectraCount = spectraCount;
        this.isUnique = isUnique;
        this.confidence = confidence;
    }

    /** default constructor */
    public DTASelectBestPeptide() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getPeptideHitId() {
        return this.peptideHitId;
    }

    public void setPeptideHitId(String peptideHitId) {
        this.peptideHitId = peptideHitId;
    }

    public String getProteinHitId() {
        return this.proteinHitId;
    }

    public void setProteinHitId(String proteinHitId) {
        this.proteinHitId = proteinHitId;
    }

    public String getSpectraCount() {
        return this.spectraCount;
    }

    public void setSpectraCount(String spectraCount) {
        this.spectraCount = spectraCount;
    }

    public Boolean getIsUnique() {
        return this.isUnique;
    }

    public void setIsUnique(Boolean isUnique) {
        this.isUnique = isUnique;
    }

    public String getConfidence() {
        return this.confidence;
    }

    public void setConfidence(String confidence) {
        this.confidence = confidence;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof DTASelectBestPeptide) ) return false;
        DTASelectBestPeptide castOther = (DTASelectBestPeptide) other;
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
