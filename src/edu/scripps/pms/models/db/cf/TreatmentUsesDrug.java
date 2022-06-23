package edu.scripps.pms.models.db.cf;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class TreatmentUsesDrug implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String drugId;

    /** nullable persistent field */
    private String treatmentId;

    /** nullable persistent field */
    private String amount;

    /** nullable persistent field */
    private String unit;

    /** nullable persistent field */
    private String note;

    /** full constructor */
    public TreatmentUsesDrug(String drugId, String treatmentId, String amount, String unit, String note) {
        this.drugId = drugId;
        this.treatmentId = treatmentId;
        this.amount = amount;
        this.unit = unit;
        this.note = note;
    }

    /** default constructor */
    public TreatmentUsesDrug() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getDrugId() {
        return this.drugId;
    }

    public void setDrugId(String drugId) {
        this.drugId = drugId;
    }

    public String getTreatmentId() {
        return this.treatmentId;
    }

    public void setTreatmentId(String treatmentId) {
        this.treatmentId = treatmentId;
    }

    public String getAmount() {
        return this.amount;
    }

    public void setAmount(String amount) {
        this.amount = amount;
    }

    public String getUnit() {
        return this.unit;
    }

    public void setUnit(String unit) {
        this.unit = unit;
    }

    public String getNote() {
        return this.note;
    }

    public void setNote(String note) {
        this.note = note;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof TreatmentUsesDrug) ) return false;
        TreatmentUsesDrug castOther = (TreatmentUsesDrug) other;
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
