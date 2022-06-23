package edu.scripps.pms.util;

import java.io.*;
import java.util.Enumeration;
import java.util.Properties;

/**
 * @author  Robin Park
 * @version $Id: AppProperties.java,v 1.3 2005/02/15 18:00:44 rpark Exp $
 */
public class AppProperties
{
    /*
    //  protected static final String
      private static Properties defaultProps = null;
      private static Properties appProps = null;
      private static String appPropsFilename = null;

      public synchronized static void initProperties(String appPropsFile) throws
          Exception
      {
          InputStream in = null;
          try {
              defaultProps = new Properties();
              Class c = AppProperties.class;
              ClassLoader cl = c.getClassLoader();
              java.net.URL sysPropURL = cl.getResource("config/sys.properties");
  //      java.net.URL sysPropURL = defaultProps.getClass().getClassLoader().getResource(
  //          Constants.SYSTEM_PROPERTY_FILE);
              in = sysPropURL.openStream();
              defaultProps.load(in);
              in.close();
          }
          catch (Throwable t) {
              t.printStackTrace();
              throw new Exception("Unable to initialize default properties!");
          }

          if (defaultProps != null) {
              appProps = new Properties(defaultProps);
          }
          else {
              appProps = new Properties();

          }
          appPropsFilename = (appPropsFile == null || appPropsFile.length() == 0) ?
              generateAppPropsFilename() : appPropsFile;

  //System.out.println("appPropsFilename = " + appPropsFilename);
          try {
              in = new FileInputStream(appPropsFilename);
              appProps.load(in);
              in.close();
          }
          catch (Throwable t) {
  //      throw new Exception ("Unable to initialize user properties!");
          }
      }

      public static Enumeration getPropertyNames()
      {
          try {
              if (appProps == null) {
                  initProperties(null);
              }
              return appProps.propertyNames();
          }
          catch (Throwable t) {
              return null;
          }
      }

      public static String getProperty(String name)
      {
          return getProperty(name, "");
      }

      public static String getProperty(String name, String defaultValue)
      {
          try {
              if (appProps == null) {
                  initProperties(null);
              }
              return appProps.getProperty(name, defaultValue);
          }
          catch (Throwable t) {
              return defaultValue;
          }
      }

      public synchronized static Object setProperty(String name, String value) throws
          Exception
      {
          if (appProps == null) {
              initProperties(null);
          }
          return appProps.setProperty(name, value);
      }

      public synchronized static void store() throws Exception
      {
  //System.out.println("Store properties to file: " + appPropsFilename);
          if (appPropsFilename == null || appProps == null) {
              return;
          }
          try {
              FileOutputStream out = new FileOutputStream(appPropsFilename);
              appProps.store(out, "-- Registration Editor Properties --");
              out.close();
          }
          catch (Throwable t) {
              throw new Exception("Unable to store application properties!");
          }
      }

      public static void store(OutputStream out, String header) throws Exception
      {
          if (appProps == null) {
              return;
          }
          try {
              appProps.store(out, header);
          }
          catch (Throwable t) {
              throw new Exception("Unable to store application properties!");
          }
      }

      private static String generateAppPropsFilename()
      {
          String user = System.getProperty("user.name");
          if (user == null) {
              user = Constants.DEFAULT_USER_NAME;
          }
          return Constants.DEFAULT_CONFIG_DIR
              + System.getProperty("file.separator")
              + user + ".properties";
      }

      private static String getDefaultUserPropsFilename()
      {
          return Constants.DEFAULT_CONFIG_DIR
              + System.getProperty("file.separator")
              + Constants.DEFAULT_USER_NAME + ".properties";

      }
*/
}
