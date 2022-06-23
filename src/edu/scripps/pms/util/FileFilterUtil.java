package edu.scripps.pms.util;

import java.io.FileFilter;
import java.io.File;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version $Id: FileFilterUtil.java,v 1.4 2005/09/20 16:54:00 rpark Exp $
 */
public class FileFilterUtil
{
    public static FileFilter getDirectoryFilter()
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return file.isDirectory();
            }
        };
    }

    public static FileFilter getFileFilter()
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return !file.isDirectory();
                // file.isFile() does not display link files
//                return file.isFile();
            }
        };
    }

    public static FileFilter getSQTFilter()
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return file.getName().endsWith("sqt");
            }
        };
    }

    public static FileFilter getExactFileFilter(final String fileName)
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return file.getName().equals(fileName);
            }
        };
    }
}
