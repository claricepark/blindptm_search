package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class SeparationType implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String separationType;

    /** full constructor */
    public SeparationType(String separationType) {
        this.separationType = separationType;
    }

    /** default constructor */
    public SeparationType() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSeparationType() {
        return this.separationType;
    }

    public void setSeparationType(String separationType) {
        this.separationType = separationType;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof SeparationType) ) return false;
        SeparationType castOther = (SeparationType) other;
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
