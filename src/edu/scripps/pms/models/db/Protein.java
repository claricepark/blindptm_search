package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Protein implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String accession;

    /** nullable persistent field */
    private String version;

    /** nullable persistent field */
    private String description;

    /** nullable persistent field */
    private String sequence;

    /** nullable persistent field */
    private String predictedPi;

    /** nullable persistent field */
    private String predictedMass;

    /** nullable persistent field */
    private String geneName;

    /** nullable persistent field */
    private String synonyms;

    /** nullable persistent field */
    private String proteinDbId;

    /** nullable persistent field */
    private edu.scripps.pms.models.db.ProteinDb proteinDb;

    /** full constructor */
    public Protein(String accession, String version, String description, String sequence, String predictedPi, String predictedMass, String geneName, String synonyms, String proteinDbId, edu.scripps.pms.models.db.ProteinDb proteinDb) {
        this.accession = accession;
        this.version = version;
        this.description = description;
        this.sequence = sequence;
        this.predictedPi = predictedPi;
        this.predictedMass = predictedMass;
        this.geneName = geneName;
        this.synonyms = synonyms;
        this.proteinDbId = proteinDbId;
        this.proteinDb = proteinDb;
    }

    /** default constructor */
    public Protein() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getAccession() {
        return this.accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getVersion() {
        return this.version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getSequence() {
        return this.sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getPredictedPi() {
        return this.predictedPi;
    }

    public void setPredictedPi(String predictedPi) {
        this.predictedPi = predictedPi;
    }

    public String getPredictedMass() {
        return this.predictedMass;
    }

    public void setPredictedMass(String predictedMass) {
        this.predictedMass = predictedMass;
    }

    public String getGeneName() {
        return this.geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getSynonyms() {
        return this.synonyms;
    }

    public void setSynonyms(String synonyms) {
        this.synonyms = synonyms;
    }

    public String getProteinDbId() {
        return this.proteinDbId;
    }

    public void setProteinDbId(String proteinDbId) {
        this.proteinDbId = proteinDbId;
    }

    public edu.scripps.pms.models.db.ProteinDb getProteinDb() {
        return this.proteinDb;
    }

    public void setProteinDb(edu.scripps.pms.models.db.ProteinDb proteinDb) {
        this.proteinDb = proteinDb;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Protein) ) return false;
        Protein castOther = (Protein) other;
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
