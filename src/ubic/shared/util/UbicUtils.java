/**
 * @(#)UbicUtils.java
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

package ubic.shared.util;

import java.io.*;
import java.util.*; // for TimeZone
import org.biojava.bio.seq.*;

/**
 * UbicUtils
 *
 * Auth Tao Xu
 * Date: November, 2003
 */

public class UbicUtils {

    /**
     * Returns the current time as a String in the format yyyy-MM-dd HH:mm:ss.
     * @return the current time as a string
     */
    public static String  time() {

        // on some JDK, the default TimeZone is wrong, we must set the TimeZone manually!!!
        //Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("EST"));

        Calendar cal = Calendar.getInstance(TimeZone.getDefault());

        String DATE_FORMAT = "yyyy-MM-dd HH:mm:ss";
        java.text.SimpleDateFormat sdf = new java.text.SimpleDateFormat(DATE_FORMAT);


        sdf.setTimeZone(TimeZone.getDefault());

        return sdf.format(cal.getTime());
    }
    
    /**
     *  Map from one strand format to another
     */
    public static final int SM_POSITIVE = 0;
    public static final int SM_NEGATIVE = 1;
    public static final int SM_UNKNOWN  = 2;
    
    public static final StrandedFeature.Strand BIOJAVA_STRAND_MAP[] = {
        StrandedFeature.POSITIVE, 
        StrandedFeature.NEGATIVE, 
        StrandedFeature.UNKNOWN
    };
    public static final int IDB_STRAND_MAP[] = {1, 2, 0};
    public static final int PEGASYS_STRAND_MAP[] = {
        StrandedFeature.POSITIVE.getValue(), 
        StrandedFeature.NEGATIVE.getValue(), 
        StrandedFeature.UNKNOWN.getValue()
    }; 
    public static final char TOKEN_STRAND_MAP[] = {
        StrandedFeature.POSITIVE.getToken(), 
        StrandedFeature.NEGATIVE.getToken(), 
        StrandedFeature.UNKNOWN.getToken()
    }; 
    
    /**
     * Functions for transforming one type of strand identifier to another. 
     * Supported strand identifiers include: the IDB Location table, the 
     * BioJava StrandedFeature.Strand object, a token (char), and the Pegasys 
     * pegasys_result table.
     * @param strand the strand field to be converted
     * @return the corresponding strand in the new format
     */        
    public static StrandedFeature.Strand idbStrand2BioJavaStrand(int strand) {
        return (strand == IDB_STRAND_MAP[SM_POSITIVE]) ? 
                      BIOJAVA_STRAND_MAP[SM_POSITIVE] : 
               (strand == IDB_STRAND_MAP[SM_NEGATIVE]) ? 
                      BIOJAVA_STRAND_MAP[SM_NEGATIVE] : 
                      BIOJAVA_STRAND_MAP[SM_UNKNOWN];
    }

    public static StrandedFeature.Strand pegasysStrand2BioJavaStrand(int strand) {
        return (strand == PEGASYS_STRAND_MAP[SM_POSITIVE]) ? 
                          BIOJAVA_STRAND_MAP[SM_POSITIVE] : 
               (strand == PEGASYS_STRAND_MAP[SM_NEGATIVE]) ? 
                          BIOJAVA_STRAND_MAP[SM_NEGATIVE] : 
                          BIOJAVA_STRAND_MAP[SM_UNKNOWN];
    }

    public static StrandedFeature.Strand tokenStrand2BioJavaStrand(char strand) {
        return (strand == TOKEN_STRAND_MAP[SM_POSITIVE]) ? 
                        BIOJAVA_STRAND_MAP[SM_POSITIVE] : 
               (strand == TOKEN_STRAND_MAP[SM_NEGATIVE]) ? 
                        BIOJAVA_STRAND_MAP[SM_NEGATIVE] : 
                        BIOJAVA_STRAND_MAP[SM_UNKNOWN];
    }

    public static int bioJavaStrand2IdbStrand(StrandedFeature.Strand strand) {
        return (strand.equals(BIOJAVA_STRAND_MAP[SM_POSITIVE])) ? 
                                  IDB_STRAND_MAP[SM_POSITIVE] : 
               (strand.equals(BIOJAVA_STRAND_MAP[SM_NEGATIVE])) ? 
                                  IDB_STRAND_MAP[SM_NEGATIVE] : 
                                  IDB_STRAND_MAP[SM_UNKNOWN];
    }

    public static int bioJavaStrand2PegasysStrand(StrandedFeature.Strand strand) {
        return strand.getValue();
    }

    public static char bioJavaStrand2TokenStrand(StrandedFeature.Strand strand) {
        return strand.getToken();
    }
    
    /**
     * Helper function to get the message and stack trace of an exception 
     * as a string (e.g. to be sent a part of a larger string).
     * @param e The exception whose information is desired.
     */
    public static String stackTrace2String(Exception e) {
        StringWriter strWriter = new StringWriter();
        e.printStackTrace(new PrintWriter(strWriter));
        return strWriter.toString();
    }
    
    /**
     * Create a temporary directory.  The java.io.File class only supports 
     * creating temporary files, not directories.
     * @param prefix A string to prepend to the directory name
     * @param suffic A string to append to the directory name
     * @param parentDir The abstract pathname of the parent directory
     * @return The newly created directory
     */
    public static File makeTempDir(String prefix, String suffix, File parentDir) {
        File dir = null;
//         while(true) {
//             // Generate a filename using the current time
//             Long time = new Long(System.currentTimeMillis());
//             String dirName = prefix + time.toString() + suffix;

//             // Make the directory
//             dir = new File(parentDir, dirName);
//             if(!dir.exists()) {
//                 dir.mkdir();
//                 break;
//             }

//             // If the name already exists, wait and try again
//             try {
//                 Thread.sleep(1000);
//             } catch (InterruptedException e) {
//                 // Ignore
//             }
//         }


	Long time = new Long(System.currentTimeMillis());
	String dirName = prefix + time.toString() + suffix;
	    
	// Make the directory
	dir = new File(parentDir, dirName);
	if(!dir.exists()) {
	    dir.mkdir();
	}
	try {
	    Thread.sleep(1000);
	}
	catch (InterruptedException e) {
	}
        return dir;
    }

    /**
     * Wrapper functions for getting/setting the working directory, 
     * in case we want to modify their behavior.
     */
    
    /**
     * Working directory setter.
     * @param workingDirectory the new working directory
     * @return void
     */
    public static void setWorkingDirectory(String workingDirectory) {
        // Changing this property doesn't affect the user directory 
        // after the program finishes
        System.setProperty("user.dir", workingDirectory);
    }
    
    /**
     * Working directory getter.
     * @return The working directory as a String
     */
    public static String getWorkingDirectory() {
        return System.getProperty("user.dir");
    }
    
    /**
     * Get the default working directory as a File object.
     * @return The working directory as a File object.
     */
    public static File getWorkingDirectoryAsFile() {
        return new File(getWorkingDirectory());
    }

    /**
     * Read a file in as a long string.
     * @param filename The name of the file to read in
     * @return The file as a String
     */
    private static String readInFile(String filename) throws IOException {
        File inputFile = new File(filename);
        int length = (int) inputFile.length(); // file must be < 2GB
        byte[] bytes = new byte[length];
        FileInputStream is = new FileInputStream(inputFile);
        is.read(bytes, 0, length);
        is.close();
        return new String(bytes);
    }

}
