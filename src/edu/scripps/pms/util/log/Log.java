package edu.scripps.pms.util.log;

//import net.nbirn.mediator.utils.Configuration;
import javax.naming.Context;
import javax.naming.InitialContext;
import javax.sql.DataSource;
import java.sql.Connection;
import java.sql.SQLException;
import java.sql.DriverManager;
//import org.jdom.Element;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.apache.log4j.Level;

//import net.nbirn.mediator.utils.EnvVarsExtractor;
import java.io.*;

import org.apache.log4j.Layout;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.FileAppender;
import java.util.Properties;
import java.util.Iterator;
import java.util.List;
import edu.scripps.pms.helper.PMSConstants;
import edu.scripps.pms.util.EnvVarsExtractor;
/**
 * A wrapper class for logging system.  Type of messages you want to pass
 * can be found in MessageType class.
 *
 * Logging requests are made by invoking one of the static methods.  A logging request is
 * enabled if its level is higher than or equal to the level of its logger.  Otherwise, the
 * request is disabled.  <p>
 * This is a basic logging rule. <br>
 * A log request of level p in a logger with level q, is enabled if  p >= q.<br>
 * Standard logging level is DEBUG < INFO < WARN < ERROR < FATAL. <p>
 *
 * Logging level and other options are defined in log4j.properties file. <p>
 *
 * Logging format is <br>
 * [DATE TIME] [Client info] [Thread number] [Package and Class name] [Method : Line number]
 * [Message Type : Detailed message]
 *
 * @author rpark
 * @version $Id: Log.java,v 1.5 2005/07/08 22:22:25 rpark Exp $
 *
 */

public class Log
{
    public static void main(String args[]) throws IOException
    {
        EnvVarsExtractor env = new EnvVarsExtractor();

        StringBuffer sb = new StringBuffer();
        sb.append(env.var(PMSConstants.APP_HOME));
        sb.append(File.separator);
        sb.append("conf");
        sb.append(File.separator);
        sb.append("log.conf");

     //   Log l = new Log(sb.toString());
        //	l.config();

//                l.debug("fdsafdsa");
    }

    private static String FQN = Log.class.getName();
    private static Logger log = Logger.getLogger(FQN);

    public static void debug(String message)
    {
        log.log(FQN, Level.DEBUG, message, null);
    }

    public static void info(String message)
    {
        log.log(FQN, Level.INFO, message, null);
    }

    public static void error(String message)
    {
        log.log(FQN, Level.ERROR, message, null);
    }

    public static void error(Exception e)
    {
	StringBuffer sb = new StringBuffer();
	
	sb.append(e.toString()).append("\n");

        StackTraceElement[] stack = e.getStackTrace();
        for(int i=0;i<stack.length;i++)
            sb.append( stack[i].toString() ).append("\n");

        error(sb.toString());
    }

    public static void fatal(String message)
    {
        log.log(FQN, Level.FATAL, message, null);
    }

    public static void warn(String message)
    {
        log.log(FQN, Level.WARN, message, null);
    }


    /*
        private static FileAppender app;
        //private String confFile;

        public Log(String confFile)
        {
            System.out.println("log start=================>>");
            PropertyConfigurator.configure(confFile);
        }

        public static void close()
        {
            app.close();
        }
    */

}
