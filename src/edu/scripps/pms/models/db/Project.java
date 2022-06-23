package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Set;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Project implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String name;

    /** nullable persistent field */
    private String description;

    /** persistent field */
    private Set mspExperiment;

    /** full constructor */
    public Project(String name, String description, Set mspExperiment) {
        this.name = name;
        this.description = description;
        this.mspExperiment = mspExperiment;
    }

    /** default constructor */
    public Project() {
    }

    /** minimal constructor */
    public Project(Set mspExperiment) {
        this.mspExperiment = mspExperiment;
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

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Set getMspExperiment() {
        return this.mspExperiment;
    }

    public void setMspExperiment(Set mspExperiment) {
        this.mspExperiment = mspExperiment;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Project) ) return false;
        Project castOther = (Project) other;
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
