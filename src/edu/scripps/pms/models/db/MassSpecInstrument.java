package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class MassSpecInstrument implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String description;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.InstrumentType instrumentTypeId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.IonSource ionSourceId;

    /** full constructor */
    public MassSpecInstrument(String description, edu.scripps.pms.models.db.InstrumentType instrumentTypeId, edu.scripps.pms.models.db.IonSource ionSourceId) {
        this.description = description;
        this.instrumentTypeId = instrumentTypeId;
        this.ionSourceId = ionSourceId;
    }

    /** default constructor */
    public MassSpecInstrument() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public edu.scripps.pms.models.db.InstrumentType getInstrumentTypeId() {
        return this.instrumentTypeId;
    }

    public void setInstrumentTypeId(edu.scripps.pms.models.db.InstrumentType instrumentTypeId) {
        this.instrumentTypeId = instrumentTypeId;
    }

    public edu.scripps.pms.models.db.IonSource getIonSourceId() {
        return this.ionSourceId;
    }

    public void setIonSourceId(edu.scripps.pms.models.db.IonSource ionSourceId) {
        this.ionSourceId = ionSourceId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof MassSpecInstrument) ) return false;
        MassSpecInstrument castOther = (MassSpecInstrument) other;
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
