package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class DbSearchParametersHasProteinDb implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String dbSearchParametersId;

    /** nullable persistent field */
    private String proteinDbId;

    /** full constructor */
    public DbSearchParametersHasProteinDb(String dbSearchParametersId, String proteinDbId) {
        this.dbSearchParametersId = dbSearchParametersId;
        this.proteinDbId = proteinDbId;
    }

    /** default constructor */
    public DbSearchParametersHasProteinDb() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getDbSearchParametersId() {
        return this.dbSearchParametersId;
    }

    public void setDbSearchParametersId(String dbSearchParametersId) {
        this.dbSearchParametersId = dbSearchParametersId;
    }

    public String getProteinDbId() {
        return this.proteinDbId;
    }

    public void setProteinDbId(String proteinDbId) {
        this.proteinDbId = proteinDbId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof DbSearchParametersHasProteinDb) ) return false;
        DbSearchParametersHasProteinDb castOther = (DbSearchParametersHasProteinDb) other;
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
