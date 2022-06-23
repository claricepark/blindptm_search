package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class ProteinHit implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String resultSelectionId;

    /** nullable persistent field */
    private String numProteins;

    /** nullable persistent field */
    private String numPeptideHit;

    /** full constructor */
    public ProteinHit(String resultSelectionId, String numProteins, String numPeptideHit) {
        this.resultSelectionId = resultSelectionId;
        this.numProteins = numProteins;
        this.numPeptideHit = numPeptideHit;
    }

    /** default constructor */
    public ProteinHit() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getResultSelectionId() {
        return this.resultSelectionId;
    }

    public void setResultSelectionId(String resultSelectionId) {
        this.resultSelectionId = resultSelectionId;
    }

    public String getNumProteins() {
        return this.numProteins;
    }

    public void setNumProteins(String numProteins) {
        this.numProteins = numProteins;
    }

    public String getNumPeptideHit() {
        return this.numPeptideHit;
    }

    public void setNumPeptideHit(String numPeptideHit) {
        this.numPeptideHit = numPeptideHit;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof ProteinHit) ) return false;
        ProteinHit castOther = (ProteinHit) other;
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
