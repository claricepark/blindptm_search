package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.ToStringBuilder;


/** @author Hibernate CodeGenerator */
public class ProteinHitHasPeptideHit implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String proteinHitId;

    /** nullable persistent field */
    private String peptideHitId;

    /** full constructor */
    public ProteinHitHasPeptideHit(String proteinHitId, String peptideHitId) {
        this.proteinHitId = proteinHitId;
        this.peptideHitId = peptideHitId;
    }

    /** default constructor */
    public ProteinHitHasPeptideHit() {
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

    public String getPeptideHitId() {
        return this.peptideHitId;
    }

    public void setPeptideHitId(String peptideHitId) {
        this.peptideHitId = peptideHitId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

}
