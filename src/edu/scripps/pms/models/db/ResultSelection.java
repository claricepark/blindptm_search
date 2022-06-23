package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class ResultSelection implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String criteria;

    /** nullable persistent field */
    private String resultSummary;

    /** nullable persistent field */
    private String numRedundantProteins;

    /** nullable persistent field */
    private String numNonredundantProteins;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.DbSearch dbSearchId;

    /** full constructor */
    public ResultSelection(String criteria, String resultSummary, String numRedundantProteins, String numNonredundantProteins, edu.scripps.pms.models.db.DbSearch dbSearchId) {
        this.criteria = criteria;
        this.resultSummary = resultSummary;
        this.numRedundantProteins = numRedundantProteins;
        this.numNonredundantProteins = numNonredundantProteins;
        this.dbSearchId = dbSearchId;
    }

    /** default constructor */
    public ResultSelection() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getCriteria() {
        return this.criteria;
    }

    public void setCriteria(String criteria) {
        this.criteria = criteria;
    }

    public String getResultSummary() {
        return this.resultSummary;
    }

    public void setResultSummary(String resultSummary) {
        this.resultSummary = resultSummary;
    }

    public String getNumRedundantProteins() {
        return this.numRedundantProteins;
    }

    public void setNumRedundantProteins(String numRedundantProteins) {
        this.numRedundantProteins = numRedundantProteins;
    }

    public String getNumNonredundantProteins() {
        return this.numNonredundantProteins;
    }

    public void setNumNonredundantProteins(String numNonredundantProteins) {
        this.numNonredundantProteins = numNonredundantProteins;
    }

    public edu.scripps.pms.models.db.DbSearch getDbSearchId() {
        return this.dbSearchId;
    }

    public void setDbSearchId(edu.scripps.pms.models.db.DbSearch dbSearchId) {
        this.dbSearchId = dbSearchId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof ResultSelection) ) return false;
        ResultSelection castOther = (ResultSelection) other;
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
