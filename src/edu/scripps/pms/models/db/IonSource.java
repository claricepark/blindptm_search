package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class IonSource implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String ionSourceType;

    /** full constructor */
    public IonSource(String ionSourceType) {
        this.ionSourceType = ionSourceType;
    }

    /** default constructor */
    public IonSource() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getIonSourceType() {
        return this.ionSourceType;
    }

    public void setIonSourceType(String ionSourceType) {
        this.ionSourceType = ionSourceType;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof IonSource) ) return false;
        IonSource castOther = (IonSource) other;
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
