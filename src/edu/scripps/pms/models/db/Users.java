package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Users implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String firstName;

    /** nullable persistent field */
    private String lastName;

    /** nullable persistent field */
    private String middleName;

    /** nullable persistent field */
    private String userName;

    /** nullable persistent field */
    private String researchGroupId;

    /** nullable persistent field */
    private String pmsUserRoleId;

    /** nullable persistent field */
    private String dataPath;

    /** full constructor */
    public Users(String firstName, String lastName, String middleName, String userName, String researchGroupId, String pmsUserRoleId, String dataPath) {
        this.firstName = firstName;
        this.lastName = lastName;
        this.middleName = middleName;
        this.userName = userName;
        this.researchGroupId = researchGroupId;
        this.pmsUserRoleId = pmsUserRoleId;
        this.dataPath = dataPath;
    }

    /** default constructor */
    public Users() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getFirstName() {
        return this.firstName;
    }

    public void setFirstName(String firstName) {
        this.firstName = firstName;
    }

    public String getLastName() {
        return this.lastName;
    }

    public void setLastName(String lastName) {
        this.lastName = lastName;
    }

    public String getMiddleName() {
        return this.middleName;
    }

    public void setMiddleName(String middleName) {
        this.middleName = middleName;
    }

    public String getUserName() {
        return this.userName;
    }

    public void setUserName(String userName) {
        this.userName = userName;
    }

    public String getResearchGroupId() {
        return this.researchGroupId;
    }

    public void setResearchGroupId(String researchGroupId) {
        this.researchGroupId = researchGroupId;
    }

    public String getPmsUserRoleId() {
        return this.pmsUserRoleId;
    }

    public void setPmsUserRoleId(String pmsUserRoleId) {
        this.pmsUserRoleId = pmsUserRoleId;
    }

    public String getDataPath() {
        return this.dataPath;
    }

    public void setDataPath(String dataPath) {
        this.dataPath = dataPath;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Users) ) return false;
        Users castOther = (Users) other;
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
