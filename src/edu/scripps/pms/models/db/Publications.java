package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Date;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Publications implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String title;

    /** nullable persistent field */
    private String path;

    /** nullable persistent field */
    private String docType;

    /** nullable persistent field */
    private String usersId;

    /** nullable persistent field */
    private Date publicationDate;

    /** nullable persistent field */
    private Date inputDate;

    /** full constructor */
    public Publications(String title, String path, String docType, String usersId, Date publicationDate, Date inputDate) {
        this.title = title;
        this.path = path;
        this.docType = docType;
        this.usersId = usersId;
        this.publicationDate = publicationDate;
        this.inputDate = inputDate;
    }

    /** default constructor */
    public Publications() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getTitle() {
        return this.title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public String getPath() {
        return this.path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getDocType() {
        return this.docType;
    }

    public void setDocType(String docType) {
        this.docType = docType;
    }

    public String getUsersId() {
        return this.usersId;
    }

    public void setUsersId(String usersId) {
        this.usersId = usersId;
    }

    public Date getPublicationDate() {
        return this.publicationDate;
    }

    public void setPublicationDate(Date publicationDate) {
        this.publicationDate = publicationDate;
    }

    public Date getInputDate() {
        return this.inputDate;
    }

    public void setInputDate(Date inputDate) {
        this.inputDate = inputDate;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Publications) ) return false;
        Publications castOther = (Publications) other;
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
