package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Date;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class ScriptsRepository implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String title;

    /** nullable persistent field */
    private String command;

    /** nullable persistent field */
    private String description;

    /** nullable persistent field */
    private Date inputDate;

    /** nullable persistent field */
    private String filename;

    /** nullable persistent field */
    private String filePath;

    /** nullable persistent field */
    private String usersId;

    /** full constructor */
    public ScriptsRepository(String title, String command, String description, Date inputDate, String filename, String filePath, String usersId) {
        this.title = title;
        this.command = command;
        this.description = description;
        this.inputDate = inputDate;
        this.filename = filename;
        this.filePath = filePath;
        this.usersId = usersId;
    }

    /** default constructor */
    public ScriptsRepository() {
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

    public String getCommand() {
        return this.command;
    }

    public void setCommand(String command) {
        this.command = command;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Date getInputDate() {
        return this.inputDate;
    }

    public void setInputDate(Date inputDate) {
        this.inputDate = inputDate;
    }

    public String getFilename() {
        return this.filename;
    }

    public void setFilename(String filename) {
        this.filename = filename;
    }

    public String getFilePath() {
        return this.filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

    public String getUsersId() {
        return this.usersId;
    }

    public void setUsersId(String usersId) {
        this.usersId = usersId;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof ScriptsRepository) ) return false;
        ScriptsRepository castOther = (ScriptsRepository) other;
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
