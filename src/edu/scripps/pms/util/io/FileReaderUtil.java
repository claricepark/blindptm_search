package edu.scripps.pms.util.io;

import java.io.File;
import java.util.Iterator;
import java.util.Arrays;

import edu.scripps.pms.util.FileFilterUtil;

/**
 * @author  Robin Park
 * @version $Id: FileReaderUtil.java,v 1.6 2006/10/25 20:33:38 rpark Exp $
 */

public class FileReaderUtil {
    public static Iterator getAllFiles(String path)
    {
        File f = new File(path);
        return Arrays.asList( f.list() ).iterator();
    }

    //Return a single SQT file only
    public static String getSQTFile(String path)
    {
        System.out.println("path==>>" + path);

        File f = new File(path.trim());

        File[] sqtFile = f.listFiles(FileFilterUtil.getSQTFilter());

        return (sqtFile.length>0)?sqtFile[0].getName():null;
    }
}
