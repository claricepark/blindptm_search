package edu.scripps.pms.util;

/**
 * @author  Robin Park
 * @version $Id: EnvVarsExtractor.java,v 1.3 2005/02/15 18:00:44 rpark Exp $
 */

import java.util.Properties;
import java.io.*;

public class EnvVarsExtractor
{

    private static Properties envVars;

    static
    {
        try {
            envVars = new Properties();

            // Find out what OS we are on.
            String osName = System.getProperty("os.name");

            // Get environment variables.
            if ("Linux".equals(osName))
                envVars.load(Runtime.getRuntime().exec("env").getInputStream());
            else
                envVars.load(
                    Runtime.getRuntime().exec("cmd.exe /C set").getInputStream()
                    );
        }
        catch (IOException ioe) {
            System.out.println("Error ==>>" + ioe);
        }
    }

    public static String var(String key)
    {
        return envVars.getProperty(key);
    }

    public static Properties envVars()
    {
        return envVars;
    }

    /*
  public static void main(String[] args) throws Exception
  {
    EnvVarsExtractor eve = new EnvVarsExtractor();
    String str = eve.var("CLASSPATH");
    System.out.println(str);
  }
 */
}
