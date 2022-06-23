/**
 * @(#)UbicConstants.java
 *
 * Copyright Notice:
 *
 * Copyright 2003 UBC Bioinformatics Centre
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

package ubic.shared.constants;

import java.util.*; // for TimeZone


/**
 * UbicConstants
 *
 * Auth Tao Xu
 * Date: November, 2003
 */

public class UbicConstants {

    /**
     * Default value for invalid integer ids
     */
    public static final int INVALIDINTID = -1;

    /**
     * Default value for invalid index, e.g., position in a String, array
     */
    public static final int INVALIDINDEX = -1;
    
    
    /**
     * Default int return value, e.g., length
     */
    public static final int INVALIDINTVALUE = -1;

    /**
     * Default value for invalid string ids
     */
    public static final String INVALIDSTRINGID = null;
    
    
    /**
     * Empty string
     */
    public static final String EMPTYSTRING = "";

    /**
     * Default value for invalid string
     */
    public static final String INVALIDSTRING = null;

    /**
     * DNA moltype
     */
    public static final int DNAMOLTYPE = 1;

    /**
     * RNA moltype
     */
    public static final int RNAMOLTYPE = 2;

    /**
     * Protein moltype
     */
    public static final int AAMOLTYPE = 3;

    /**
     * Nucleic acid moltype
     */
    public static final int NAMOLTYPE = 4;


    public static final int ASCENDINGORDER = 1;

    public static final int DESCENDINGORDER = 2;

    public static final int STATISTICBASED = 1;
    
    public static final int SCOREBASED = 2;

    public static final double JDBCFLOATTHRESHOLD = 1.0E-99;;

}
