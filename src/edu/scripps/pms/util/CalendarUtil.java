package edu.scripps.pms.util;

import java.util.Calendar;
import java.util.*;
import java.sql.Date;
import java.text.SimpleDateFormat;

/**
 * @author Robin Park
 * @version $Id: CalendarUtil.java,v 1.5 2006/09/27 21:50:14 rpark Exp $
 */
public class CalendarUtil
{
    private Calendar cal;


    public CalendarUtil(int year, int month, int date)
    {
        init(year, month, date);
    }

    /**
     * CalendarUtil
     *
     * @param dateString String
     * Format of param is mm/dd/yyyy
     */
    public CalendarUtil(String dateString)
    {
        String[] str = dateString.split("/");

        //Fix me later.
        if( str.length != 3 )
            System.out.println("String input format is incorrect");

        //Month starts from '0', but date starts from 1.  interesting..
        init(Integer.parseInt(str[2]), Integer.parseInt(str[0])-1, Integer.parseInt(str[1]));
    }

    public void init(int year, int month, int date)
    {
        cal = Calendar.getInstance();
        cal.set(year, month, date);
    }

    public long getDateAsLong()
    {
        return cal.getTime().getTime();
    }

    public Date getDate()
    {
        return new Date( getDateAsLong() );
    }

    public static String getMediumFormat(Date date)
    {
        return TimeUtils.getMediumFormat().format(date);
    }

    public static String getDateTimeFormat(Date date)
    {
        return TimeUtils.getDateTimeFormat().format(date);
    }

    public static String getCurrentDateTime()
    {
        //sample format
        //"dd MMMMM yyyy"
        //"yyyyMMdd"
        //"dd.MM.yy"
        //"MM/dd/yy"
        //"yyyy.MM.dd G 'at' hh:mm:ss z"
        //"EEE, MMM d, ''yy"
        //"h:mm a"
        //"H:mm:ss:SSS"
        //"K:mm a,z"
        //"yyyy.MMMMM.dd GGG hh:mm aaa"

        return getCurrentDateTime("yyyyMMdd_hhmmssa");  //year month date hour min second AM/PM
    }

    public static String getCurrentDateTime(String format)
    {
        java.util.Date today = new java.util.Date();
        SimpleDateFormat formatter = new SimpleDateFormat(format);

        return formatter.format(today);
    }

    public static List getMonthList()
    {
	List l = new Vector();
	for(int i=1;i<=12;i++)
	    l.add(i);

	return l;
    }	

}
