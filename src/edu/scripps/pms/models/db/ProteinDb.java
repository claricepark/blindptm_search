package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.sql.Date;
import java.sql.Timestamp;

import edu.scripps.pms.util.CalendarUtil;

import org.apache.commons.lang.builder.ToStringBuilder;


/** @author Hibernate CodeGenerator */
public class ProteinDb implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String organismId;

    /** nullable persistent field */
    private String fileName;

    /** nullable persistent field */
    private Date releaseDate;

    /** nullable persistent field */
    private String path;

    /** nullable persistent field */
    private String description;

    /** nullable persistent field */
    private String userId;

    /** nullable persistent field */
    private Timestamp inputDate;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.DbSource dbSource;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.Organism organism;

    /** full constructor */
    public ProteinDb(String organismId, String fileName, Date releaseDate, String path, String description, String userId, Timestamp inputDate, edu.scripps.pms.models.db.DbSource dbSource, edu.scripps.pms.models.db.Organism organism) {
        this.organismId = organismId;
        this.fileName = fileName;
        this.releaseDate = releaseDate;
        this.path = path;
        this.description = description;
        this.userId = userId;
        this.inputDate = inputDate;
        this.dbSource = dbSource;
        this.organism = organism;
    }

    /** default constructor */
    public ProteinDb() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getOrganismId() {
        return this.organismId;
    }

    public void setOrganismId(String organismId) {
        this.organismId = organismId;
    }

    public String getFileName() {
        return this.fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    public Date getReleaseDate() {
        return this.releaseDate;
    }

    public void setReleaseDate(Date releaseDate) {
        this.releaseDate = releaseDate;
    }

    public String getPath() {
        return this.path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getUserId() {
        return this.userId;
    }

    public void setUserId(String userId) {
        this.userId = userId;
    }

    public edu.scripps.pms.models.db.DbSource getDbSource() {
        return this.dbSource;
    }

    public void setDbSource(edu.scripps.pms.models.db.DbSource dbSource) {
        this.dbSource = dbSource;
    }

    public edu.scripps.pms.models.db.Organism getOrganism() {
        return this.organism;
    }

    public void setOrganism(edu.scripps.pms.models.db.Organism organism) {
        this.organism = organism;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public Timestamp getInputDate() {
        return this.inputDate;
    }

    public String getInputDateTimeFormat()
    {
        return CalendarUtil.getDateTimeFormat(new Date(this.inputDate.getTime()));
    }

    public String getDbDateMediumFormat() {
        return CalendarUtil.getMediumFormat(this.releaseDate);
    }

    public void setInputDate(Timestamp inputDate) {
        this.inputDate = inputDate;
    }

}
