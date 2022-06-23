package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Set;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class PeptideHit implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String msmsfractionId;

    /** nullable persistent field */
    private String scoreTypeId;

    /** nullable persistent field */
    private String dbSearchId;

    /** nullable persistent field */
    private String sequence;

    /** nullable persistent field */
    private String predictedMass;

    /** nullable persistent field */
    private String primaryScore;

    /** nullable persistent field */
    private String primaryRank;

    /** nullable persistent field */
    private String deltaCn;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.Msmsfraction msmsFraction;

    /** persistent field */
    private Set ptmList;

    /** full constructor */
    public PeptideHit(String msmsfractionId, String scoreTypeId, String dbSearchId, String sequence, String predictedMass, String primaryScore, String primaryRank, String deltaCn, edu.scripps.pms.models.db.Msmsfraction msmsFraction, Set ptmList) {
        this.msmsfractionId = msmsfractionId;
        this.scoreTypeId = scoreTypeId;
        this.dbSearchId = dbSearchId;
        this.sequence = sequence;
        this.predictedMass = predictedMass;
        this.primaryScore = primaryScore;
        this.primaryRank = primaryRank;
        this.deltaCn = deltaCn;
        this.msmsFraction = msmsFraction;
        this.ptmList = ptmList;
    }

    /** default constructor */
    public PeptideHit() {
    }

    /** minimal constructor */
    public PeptideHit(Set ptmList) {
        this.ptmList = ptmList;
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getMsmsfractionId() {
        return this.msmsfractionId;
    }

    public void setMsmsfractionId(String msmsfractionId) {
        this.msmsfractionId = msmsfractionId;
    }

    public String getScoreTypeId() {
        return this.scoreTypeId;
    }

    public void setScoreTypeId(String scoreTypeId) {
        this.scoreTypeId = scoreTypeId;
    }

    public String getDbSearchId() {
        return this.dbSearchId;
    }

    public void setDbSearchId(String dbSearchId) {
        this.dbSearchId = dbSearchId;
    }

    public String getSequence() {
        return this.sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getPredictedMass() {
        return this.predictedMass;
    }

    public void setPredictedMass(String predictedMass) {
        this.predictedMass = predictedMass;
    }

    public String getPrimaryScore() {
        return this.primaryScore;
    }

    public void setPrimaryScore(String primaryScore) {
        this.primaryScore = primaryScore;
    }

    public String getPrimaryRank() {
        return this.primaryRank;
    }

    public void setPrimaryRank(String primaryRank) {
        this.primaryRank = primaryRank;
    }

    public String getDeltaCn() {
        return this.deltaCn;
    }

    public void setDeltaCn(String deltaCn) {
        this.deltaCn = deltaCn;
    }

    public edu.scripps.pms.models.db.Msmsfraction getMsmsFraction() {
        return this.msmsFraction;
    }

    public void setMsmsFraction(edu.scripps.pms.models.db.Msmsfraction msmsFraction) {
        this.msmsFraction = msmsFraction;
    }

    public Set getPtmList() {
        return this.ptmList;
    }

    public void setPtmList(Set ptmList) {
        this.ptmList = ptmList;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof PeptideHit) ) return false;
        PeptideHit castOther = (PeptideHit) other;
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
