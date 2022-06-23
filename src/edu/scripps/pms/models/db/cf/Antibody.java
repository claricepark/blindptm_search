package edu.scripps.pms.models.db.cf;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Antibody implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String name;

    /** nullable persistent field */
    private String source;

    /** nullable persistent field */
    private String antigenName;

    /** nullable persistent field */
    private String description;

    /** full constructor */
    public Antibody(String name, String source, String antigenName, String description) {
        this.name = name;
        this.source = source;
        this.antigenName = antigenName;
        this.description = description;
    }

    /** default constructor */
    public Antibody() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getName() {
        return this.name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getSource() {
        return this.source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public String getAntigenName() {
        return this.antigenName;
    }

    public void setAntigenName(String antigenName) {
        this.antigenName = antigenName;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Antibody) ) return false;
        Antibody castOther = (Antibody) other;
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
