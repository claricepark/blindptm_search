/*
 * @(#)GoAssociation.java
 *
 * Copyright Notice:
 *
 * Copyright 2004 UBC Bioinformatics Centre
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 *        Free Software Foundation, Inc.,
 *        59 Temple Place, Suite 330, Boston, MA
 *        02111-1307  USA
 * or visit http://www.gnu.org/copyleft/gpl.html
 *
 */

/**
 * @file GoAssociation.java
 * This is the source file for ubic.atlas.locuslink.db.GoAssociation
 *
 * @author Tao Xu
 * @date $Date
 */


package edu.scripps.pms.go;


import java.util.*;





/**
 * @class GoAssociation
 * This class is a subclass of MysqlDb. It provides an interface to
 * access GO database  
 *
 */
public class GoAssociation {
    private String line;
    private String [] arr;
    public GoAssociation(String associationLine) {
        line = associationLine;
        arr = line.split("\t");        
      
    }
    public String getDb() {
        return arr[0];
    }
    public String getDbObjectId() {
        return arr[1];
    }
    public String getDbObjectSymbol() {
        return arr[2];
    }
    public String getQualifier() {
        return arr[3];
    }
    public String getGoId() {
        return arr[4];
    }
    public String getDbObjectSynonym() {
        return arr[10];
    }
    public String getTaxon() {
        return arr[12];
    }
}





