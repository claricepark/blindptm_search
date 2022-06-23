package edu.scripps.pms.models;

import java.util.Date;

public class DbSearchParameters {
    private String program;
    private String db;
    private Date dbDate;
    private String parametersFile;
    private String parameters;
    private float peptideMassTolerance;
    private float fragmentIonTolerance;
    private int maxMissedCleavages;
    private int massTypeParent;
    private int massTypeFragment;
    private int enzymeNumber;

    public void setEnzymeNumber(int enzymeNumber) {
        this.enzymeNumber = enzymeNumber;
    }

    public void setMassTypeFragment(int massTypeFragment) {
        this.massTypeFragment = massTypeFragment;
    }

    public void setMassTypeParent(int massTypeParent) {
        this.massTypeParent = massTypeParent;
    }

    public void setMaxMissedCleavages(int maxMissedCleavages) {
        this.maxMissedCleavages = maxMissedCleavages;
    }

    public void setFragmentIonTolerance(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    public void setPeptideMassTolerance(float peptideMassTolerance) {
        this.peptideMassTolerance = peptideMassTolerance;
    }

    public void setParameters(String parameters) {
        this.parameters = parameters;
    }

    public void setParametersFile(String parametersFile) {
        this.parametersFile = parametersFile;
    }

    public void setDb(String db) {
        this.db = db;
    }

    public void setProgram(String program) {
        this.program = program;
    }

    public String getProgram() {
        return program;
    }

    public String getDb() {
        return db;
    }

    public String getParametersFile() {
        return parametersFile;
    }

    public String getParameters() {
        return parameters;
    }

    public float getPeptideMassTolerance() {
        return peptideMassTolerance;
    }

    public float getFragmentIonTolerance() {
        return fragmentIonTolerance;
    }

    public int getMaxMissedCleavages() {
        return maxMissedCleavages;
    }

    public int getMassTypeParent() {
        return massTypeParent;
    }

    public int getMassTypeFragment() {
        return massTypeFragment;
    }

    public int getEnzymeNumber() {
        return enzymeNumber;
    }

}
