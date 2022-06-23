package edu.scripps.pms.models.db.cf;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class IpProtocol implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String name;

    /** nullable persistent field */
    private String antibodyId;

    /** nullable persistent field */
    private String protocol;

    /** nullable persistent field */
    private String description;

    /** full constructor */
    public IpProtocol(String name, String antibodyId, String protocol, String description) {
        this.name = name;
        this.antibodyId = antibodyId;
        this.protocol = protocol;
        this.description = description;
    }

    /** default constructor */
    public IpProtocol() {
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

    public String getAntibodyId() {
        return this.antibodyId;
    }

    public void setAntibodyId(String antibodyId) {
        this.antibodyId = antibodyId;
    }

    public String getProtocol() {
        return this.protocol;
    }

    public void setProtocol(String protocol) {
        this.protocol = protocol;
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
        if ( !(other instanceof IpProtocol) ) return false;
        IpProtocol castOther = (IpProtocol) other;
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
