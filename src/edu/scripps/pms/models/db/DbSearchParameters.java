package edu.scripps.pms.models.db;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class DbSearchParameters implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String parameters;

    /** nullable persistent field */
    private String proteinDbId;

    /** nullable persistent field */
    private String programId;

    /** nullable persistent field */
    private String parametersFile;

    /** nullable persistent field */
    private String peptideMassTolerance;

    /** nullable persistent field */
    private String fragmentIonTolerance;

    /** nullable persistent field */
    private String maxMissedCleavages;

    /** nullable persistent field */
    private String massTypeParent;

    /** nullable persistent field */
    private String massTypeFragment;

    /** nullable persistent field */
    private String enzymeNumber;

    /** full constructor */
    public DbSearchParameters(String parameters, String proteinDbId, String programId, String parametersFile, String peptideMassTolerance, String fragmentIonTolerance, String maxMissedCleavages, String massTypeParent, String massTypeFragment, String enzymeNumber) {
        this.parameters = parameters;
        this.proteinDbId = proteinDbId;
        this.programId = programId;
        this.parametersFile = parametersFile;
        this.peptideMassTolerance = peptideMassTolerance;
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.maxMissedCleavages = maxMissedCleavages;
        this.massTypeParent = massTypeParent;
        this.massTypeFragment = massTypeFragment;
        this.enzymeNumber = enzymeNumber;
    }

    /** default constructor */
    public DbSearchParameters() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getParameters() {
        return this.parameters;
    }

    public void setParameters(String parameters) {
        this.parameters = parameters;
    }

    public String getProteinDbId() {
        return this.proteinDbId;
    }

    public void setProteinDbId(String proteinDbId) {
        this.proteinDbId = proteinDbId;
    }

    public String getProgramId() {
        return this.programId;
    }

    public void setProgramId(String programId) {
        this.programId = programId;
    }

    public String getParametersFile() {
        return this.parametersFile;
    }

    public void setParametersFile(String parametersFile) {
        this.parametersFile = parametersFile;
    }

    public String getPeptideMassTolerance() {
        return this.peptideMassTolerance;
    }

    public void setPeptideMassTolerance(String peptideMassTolerance) {
        this.peptideMassTolerance = peptideMassTolerance;
    }

    public String getFragmentIonTolerance() {
        return this.fragmentIonTolerance;
    }

    public void setFragmentIonTolerance(String fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    public String getMaxMissedCleavages() {
        return this.maxMissedCleavages;
    }

    public void setMaxMissedCleavages(String maxMissedCleavages) {
        this.maxMissedCleavages = maxMissedCleavages;
    }

    public String getMassTypeParent() {
        return this.massTypeParent;
    }

    public void setMassTypeParent(String massTypeParent) {
        this.massTypeParent = massTypeParent;
    }

    public String getMassTypeFragment() {
        return this.massTypeFragment;
    }

    public void setMassTypeFragment(String massTypeFragment) {
        this.massTypeFragment = massTypeFragment;
    }

    public String getEnzymeNumber() {
        return this.enzymeNumber;
    }

    public void setEnzymeNumber(String enzymeNumber) {
        this.enzymeNumber = enzymeNumber;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof DbSearchParameters) ) return false;
        DbSearchParameters castOther = (DbSearchParameters) other;
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
