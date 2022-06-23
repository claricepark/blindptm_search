package edu.scripps.pms.models.db;

import java.io.Serializable;
import java.util.Set;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class PeakList implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String dataSourceId;

    /** nullable persistent field */
    private String listType;

    /** nullable persistent field */
    private String loscan;

    /** nullable persistent field */
    private String spectrum;

    /** nullable persistent field */
    private String mtoz;

    /** persistent field */
    private Set msmsfractions;

    /** full constructor */
    public PeakList(String dataSourceId, String listType, String loscan, String spectrum, String mtoz, Set msmsfractions) {
        this.dataSourceId = dataSourceId;
        this.listType = listType;
        this.loscan = loscan;
        this.spectrum = spectrum;
        this.mtoz = mtoz;
        this.msmsfractions = msmsfractions;
    }

    /** default constructor */
    public PeakList() {
    }

    /** minimal constructor */
    public PeakList(Set msmsfractions) {
        this.msmsfractions = msmsfractions;
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getDataSourceId() {
        return this.dataSourceId;
    }

    public void setDataSourceId(String dataSourceId) {
        this.dataSourceId = dataSourceId;
    }

    public String getListType() {
        return this.listType;
    }

    public void setListType(String listType) {
        this.listType = listType;
    }

    public String getLoscan() {
        return this.loscan;
    }

    public void setLoscan(String loscan) {
        this.loscan = loscan;
    }

    public String getSpectrum() {
        return this.spectrum;
    }

    public void setSpectrum(String spectrum) {
        this.spectrum = spectrum;
    }

    public String getMtoz() {
        return this.mtoz;
    }

    public void setMtoz(String mtoz) {
        this.mtoz = mtoz;
    }

    public Set getMsmsfractions() {
        return this.msmsfractions;
    }

    public void setMsmsfractions(Set msmsfractions) {
        this.msmsfractions = msmsfractions;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof PeakList) ) return false;
        PeakList castOther = (PeakList) other;
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
