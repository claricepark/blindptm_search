package edu.scripps.pms.util;

/**
 * @author  Robin Park
 * @version $Id: StringUtil.java,v 1.8 2007/06/19 18:21:01 rpark Exp $
 */
public class StringUtil
{
    public static String splitWithReturn(String str, int eachSize, String deliminator)
    {
        int length = str.length();
        StringBuffer sb = new StringBuffer();

        for (int i=0; i<length; i+=eachSize)
        {

            if (i + eachSize > length) {
                sb.append( str.substring(i, length) );
                sb.append(deliminator);
            }
            else {
                sb.append( str.substring(i, i+eachSize) );
                sb.append(deliminator);
            }
        }

        return sb.toString();
    }

    /* size : the location to cut of the string
	Originally this method is needed for the web application to display
	lengthy def line.  If it is too long, display just part of it with "..."
    */
    public static String trimString(String str, int size)
    {
	int length = str.length();

	boolean isOversize=false;
	if(length<=size)
		isOversize=true;

	str = str.substring(0, isOversize?length:size);
	str += isOversize?"":"...";
	return str;

    }

    /*
	remove the string before the first space.
	Originally this method is needed for the web application to display
	lengthy def line.  If it is too long, display just part of it with "..."
    */
    public static String trimAccession(String str, int size)
    {
	str = str.substring(str.indexOf(" ")+1);

	return trimString(str, size);
    }

    public static String getRandomString(int length) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < length; i++) {
            sb.append((char) (65 + (int) (Math.random() * (90 - 65))));

        }
	sb.append(System.currentTimeMillis());

        return sb.toString();
    }

    public static void main(String args[])
    {
	System.out.println( StringUtil.trimAccession("fdsafdsa99fdsaffdsafdsaffdsafdsafsdafdsdfffffffffffffffddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddfsdafds", 70) );

    }
}
