package edu.scripps.pms.util;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version $Id
 */
public class ByteArrayConverter {
    public static final int getInt(byte [] b, int index) {
        int i = 0;
        i = b[index];
        i <<= 8;
        i |= b[index + 1] & 0xFF;
        i <<= 8;
        i |= b[index + 2] & 0xFF;
        i <<= 8;
        i |= b[index + 3] & 0xFF;
        return i;
    }
    public static final int getUnsignedShort(byte [] b, int index) {
        int i = 0;
        i |= b[index] & 0xFF;
        i <<= 8;
        i |= b[index+1] & 0xFF;
        return i;
    }
    /**
     * Converts a 4 byte array of unsigned bytes to an long
     * @param b an array of 4 unsigned bytes
     * @return a long representing the unsigned int
     */
    public static final long unsignedIntToLong(byte[] b) {
        long l = 0;
        l |= b[0] & 0xFF;
        l <<= 8;
        l |= b[1] & 0xFF;
        l <<= 8;
        l |= b[2] & 0xFF;
        l <<= 8;
        l |= b[3] & 0xFF;
        return l;
    }
    
    /**
     * Converts a two byte array to an integer
     * @param b a byte array of length 2
     * @return an int representing the unsigned short
     */
    public static final int unsignedShortToInt(byte[] b) {
        int i = 0;
        i |= b[0] & 0xFF;
        i <<= 8;
        i |= b[1] & 0xFF;
        return i;
    }

}
