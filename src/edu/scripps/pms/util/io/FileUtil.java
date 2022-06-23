
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
 * @version $Id: FileUtil.java,v 1.4 2008/07/15 20:34:36 taoxu Exp $
 *
 */

package edu.scripps.pms.util.io;

import java.io.*;

import java.nio.channels.FileChannel;
import java.nio.MappedByteBuffer;

public class FileUtil
{
    private FileUtil()
    {
    }

    public static void copy(String in, String out, boolean isAppending) throws IOException {
        copy(new File(in), new File(out), isAppending);
    }

    public static void copy(File in, File out, boolean isAppending) throws IOException {
        FileInputStream fis  = new FileInputStream(in);
        FileOutputStream fos = new FileOutputStream(out, isAppending);
        byte[] buf = new byte[4096];
        int i = 0;

        while((i=fis.read(buf))!=-1) {
          fos.write(buf, 0, i);
        }

        fis.close();
        fos.close();
    }


    //This method is using new lib from jdk and supposed to be faster then old lib.
    //But before using it, make sure it is working correctly.
    //Currently we are not using this method.
    /*
    public static void copy(FileInputStream source, FileOutputStream dest) throws IOException {
         FileChannel in = null, out = null;
         try
         {
              in = source.getChannel();
              out = dest.getChannel();

              long size = in.size();
              MappedByteBuffer buf = in.map(FileChannel.MapMode.READ_ONLY, 0, size);

              out.write(buf);

         }
         catch(IOException e)
         {
             throw new IOException(e.toString());
         }
         finally {
              if (in != null)          in.close();
              if (out != null)     out.close();
         }
    }
    */
}
