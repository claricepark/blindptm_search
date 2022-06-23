package edu.scripps.pms.pathway;

import java.util.*;

public class PathwayModel implements Comparable
{
    private String name;
    private Set set = new HashSet();
    private double twoPValue;
    private double rightPValue;

    public String getName() {
	return name;
    }

    public void setName(String name) {
	this.name = name;
    }

    public Set getSet() {
	return set;
    }

    public void setSet(Set set) {
	this.set = set;
    }

    public double getTwoPValue() {
	return twoPValue;
    }

    public void setTwoPValue(double twoPValue) {
	this.twoPValue = twoPValue;
    }

    public double getRightPValue() {
	return rightPValue;
    }

    public void setRightPValue(double rightPValue) {
	this.rightPValue = rightPValue;
    }

    public int compareTo(Object obj)
    {
	//double d = (this.getTwoPValue() - ((PathwayModel)obj).getTwoPValue());
	double d = (this.getRightPValue() - ((PathwayModel)obj).getRightPValue());

	if(d>0)
	    return 1;
	else if(d==0)
	    return 0;
	else
	    return -1;
	
    }
}
