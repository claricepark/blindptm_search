package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Msmsfraction implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String peakListId;

    /** nullable persistent field */
    private String mtoz;

    /** nullable persistent field */
    private String chargeState;

    /** full constructor */
    public Msmsfraction(String peakListId, String mtoz, String chargeState) {
        this.peakListId = peakListId;
        this.mtoz = mtoz;
        this.chargeState = chargeState;
    }

    /** default constructor */
    public Msmsfraction() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getPeakListId() {
        return this.peakListId;
    }

    public void setPeakListId(String peakListId) {
        this.peakListId = peakListId;
    }

    public String getMtoz() {
        return this.mtoz;
    }

    public void setMtoz(String mtoz) {
        this.mtoz = mtoz;
    }

    public String getChargeState() {
        return this.chargeState;
    }

    public void setChargeState(String chargeState) {
        this.chargeState = chargeState;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Msmsfraction) ) return false;
        Msmsfraction castOther = (Msmsfraction) other;
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
