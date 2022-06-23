package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class ProteinHitIdentifiedProtein implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String proteinHitId;

    /** nullable persistent field */
    private String sequenceCoverage;

    /** nullable persistent field */
    private String evaluationState;

    /** nullable persistent field */
    private String confidence;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.Protein proteinId;

    /** full constructor */
    public ProteinHitIdentifiedProtein(String proteinHitId, String sequenceCoverage, String evaluationState, String confidence, edu.scripps.pms.models.db.Protein proteinId) {
        this.proteinHitId = proteinHitId;
        this.sequenceCoverage = sequenceCoverage;
        this.evaluationState = evaluationState;
        this.confidence = confidence;
        this.proteinId = proteinId;
    }

    /** default constructor */
    public ProteinHitIdentifiedProtein() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getProteinHitId() {
        return this.proteinHitId;
    }

    public void setProteinHitId(String proteinHitId) {
        this.proteinHitId = proteinHitId;
    }

    public String getSequenceCoverage() {
        return this.sequenceCoverage;
    }

    public void setSequenceCoverage(String sequenceCoverage) {
        this.sequenceCoverage = sequenceCoverage;
    }

    public String getEvaluationState() {
        return this.evaluationState;
    }

    public void setEvaluationState(String evaluationState) {
        this.evaluationState = evaluationState;
    }

    public String getConfidence() {
        return this.confidence;
    }

    public void setConfidence(String confidence) {
        this.confidence = confidence;
    }

    public edu.scripps.pms.models.db.Protein getProteinId() {
        return this.proteinId;
    }

    public void setProteinId(edu.scripps.pms.models.db.Protein proteinId) {
        this.proteinId = proteinId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof ProteinHitIdentifiedProtein) ) return false;
        ProteinHitIdentifiedProtein castOther = (ProteinHitIdentifiedProtein) other;
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
