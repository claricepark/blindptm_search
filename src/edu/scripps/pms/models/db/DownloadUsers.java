package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Date;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class DownloadUsers implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String firstName;

    /** nullable persistent field */
    private String middleName;

    /** nullable persistent field */
    private String lastName;

    /** nullable persistent field */
    private String organization;

    /** nullable persistent field */
    private String email;

    /** nullable persistent field */
    private Date inputDate;

    /** nullable persistent field */
    private String softwareId;

    /** nullable persistent field */
    private String downloadKey;

    /** nullable persistent field */
    private String ipAddress;

    /** full constructor */
    public DownloadUsers(String firstName, String middleName, String lastName, String organization, String email, Date inputDate, String softwareId, String downloadKey, String ipAddress) {
        this.firstName = firstName;
        this.middleName = middleName;
        this.lastName = lastName;
        this.organization = organization;
        this.email = email;
        this.inputDate = inputDate;
        this.softwareId = softwareId;
        this.downloadKey = downloadKey;
        this.ipAddress = ipAddress;
    }

    /** default constructor */
    public DownloadUsers() {
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

    public String getMiddleName() {
        return this.middleName;
    }

    public void setMiddleName(String middleName) {
        this.middleName = middleName;
    }

    public String getLastName() {
        return this.lastName;
    }

    public void setLastName(String lastName) {
        this.lastName = lastName;
    }

    public String getOrganization() {
        return this.organization;
    }

    public void setOrganization(String organization) {
        this.organization = organization;
    }

    public String getEmail() {
        return this.email;
    }

    public void setEmail(String email) {
        this.email = email;
    }

    public Date getInputDate() {
        return this.inputDate;
    }

    public void setInputDate(Date inputDate) {
        this.inputDate = inputDate;
    }

    public String getSoftwareId() {
        return this.softwareId;
    }

    public void setSoftwareId(String softwareId) {
        this.softwareId = softwareId;
    }

    public String getDownloadKey() {
        return this.downloadKey;
    }

    public void setDownloadKey(String downloadKey) {
        this.downloadKey = downloadKey;
    }

    public String getIpAddress() {
        return this.ipAddress;
    }

    public void setIpAddress(String ipAddress) {
        this.ipAddress = ipAddress;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof DownloadUsers) ) return false;
        DownloadUsers castOther = (DownloadUsers) other;
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
