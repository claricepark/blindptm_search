package edu.scripps.pms.util.sqt;

import java.io.BufferedReader;
import java.util.*;
import java.io.*;
//import edu.scripps.pms.util.sqt.SQTPeptide;
//import edu.scripps.pms.util.sqt.MLine;

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
 * @version $Id: SQTValidator.java,v 1.15 2006/05/01 17:15:36 rpark Exp $
 */
public class SQTValidator
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;
    private int lineNumber=1;
    private StringBuffer header = new StringBuffer();

    private static boolean fixit=false;
    private static boolean backup=false;
    private static boolean fixSingleFile=false;
    private static String singleFile="";

    public static void printUsage()
    {
            System.out.println("Usage: sqt_validator [OPTION]");
            System.out.println("-f : remove corrupted spectra from sqt files.");
            System.out.println("-s : specify file names to validate.");
            System.out.println("Example: sqt_validator -s test.sqt");
            System.out.println("Example: sqt_validator -s test.sqt -f");
    }

    public static void main(String args[]) throws IOException, Exception
    {
        String eachFile;
        if( args.length>0 && "--help".equals(args[0]) )
        {
	    printUsage();
            System.exit(0);
        }

	for(int i=0;i<args.length;i++)
	{
	    if( "-s".equals(args[i]) )
	    {
		fixSingleFile=true;
	
		if(args.length<=(i+1) || null == args[i+1] || "".equals(args[i+1]))
		{
		    System.out.println("File name to fix required");
		    printUsage();
		    System.exit(0);
		}

		singleFile=args[i+1];		
	    } else if("-f".equals(args[i]))
	    {
		fixit=true;
	    } else if("-b".equals(args[i]))
	    {
		backup=true;
	    }
	}

	if(fixSingleFile)
	{

	    if(!fixit)
	    {
		    System.out.println("Validating and fixing file " + singleFile + "...");

		    //SQTValidator reader = new SQTValidator("/data/2/taoxu/testlong/" + eachFile, args);
		    SQTValidator reader = new SQTValidator(singleFile, args);
		    reader.prevalidate();  //check if there are non H, M, L lines
		    reader.validate();
		    reader.close();
		    System.out.println("");
	    }
	    else
	    {
		    System.out.println("Validating file without fixing " + singleFile + "...");

		    //SQTValidator reader = new SQTValidator("/data/2/taoxu/testlong/" + eachFile, args);
		    SQTValidator reader = new SQTValidator(singleFile, args);
		    reader.prevalidate();  //check if there are non H, M, L lines
		    reader.validate();
		    reader.close();
		    System.out.println("");
	    }

	    System.out.println("Validating sqt files finished");
	    System.exit(0);	    
	}
	
	String path=".";
//	String path="/data/1/emilyc/Lindsey_exp/frac6/parc/";
	
        //for( Iterator fileItr = getAllFiles("/data/1/rpark/modification_cruse/pH35/MS3/ST18M48"); fileItr.hasNext(); )
        for( Iterator fileItr = getAllFiles(path); fileItr.hasNext(); )
        {
            eachFile = fileItr.next().toString();

            if(eachFile.endsWith("sqt"))
            {
                System.out.println("Validating file " + eachFile + "...");

         	//SQTValidator reader = new SQTValidator("/data/2/taoxu/testlong/" + eachFile, args);
         	SQTValidator reader = new SQTValidator(eachFile, args);
		reader.prevalidate();  //check if there are non H, M, L lines
         	reader.validate();
	 	reader.close();
                System.out.println("");

  //              System.out.println(reader.getHeader());
            }
        }

        System.out.println("Validating sqt files finished");

    }

    public void checkSinglefileFix(boolean fixit, String singleFile)
    {
    }

    public SQTValidator(String fileName, String[] args) throws IOException, Exception
    {
        this.fileName = fileName;

/*
        if(args.length==1)
        {
            if("-f".equals(args[0]) )
                fixit = true;
        }


        if(args.length==2)
        {
            if( "-f".equals(args[0]) || "-f".equals(args[1]) )
                fixit = true;

            if( "-b".equals(args[0]) || "-b".equals(args[1]) )
                backup = true;
        }


*/

    }

    private void prevalidate() throws IOException, Exception
    {
        br = new BufferedReader(new FileReader(fileName));

	boolean isError=false;

	lastLine = br.readLine();
	int index=1;

	StringBuffer fixedBuffer = new StringBuffer();
        while( lastLine.startsWith("H\t") || lastLine.startsWith(" ")) {
            fixedBuffer.append(lastLine).append("\n");
            lastLine = br.readLine(); index++;
        } //read head line

        while( lastLine.equals("") ) { lastLine = br.readLine(); index++; fixedBuffer.append("\n").append(lastLine); }//read empty line
       //Check S line

	while( (lastLine=br.readLine())!= null )
	{
	    index++;
		    if(lastLine.equals("")) continue;

	    if( !lastLine.startsWith("S\t") && !lastLine.startsWith("M\t") && !lastLine.startsWith("L\t") )
	    {
		isError=true;
		System.out.println("Error=====>> Unexpected line : Line  + # " + index + "\t" + lastLine);
	    } else {
		fixedBuffer.append("\n").append(lastLine);
	    }
	}

        //backup the sqt file
        if(!fixit)
            return;

        if(isError)
        {
            StringBuffer command = new StringBuffer();

            if(backup)
            {
                command.append("mv ");
                command.append(fileName).append(" ");
                command.append(fileName).append(".bak");
            } else
            {
                command.append("rm ");
                command.append(fileName);
            }


            Process p = Runtime.getRuntime().exec(command.toString());
            p.waitFor();

            // Connect print stream to the output stream
            PrintStream out = new PrintStream(new FileOutputStream(fileName));

            out.println(fixedBuffer.toString());

            if(backup)
            {
                System.out.print(fileName + " was fixed.  Original file was saved as ");
                System.out.println(fileName + ".bak");
            } else
            {
                System.out.print(fileName + " was fixed.");
            }

            out.close();
        }

    }


    private void validate() throws IOException, Exception
    {
        //br = new BufferedReader(new FileReader("/data/1/rpark/modification_cruse/pH35/MS3/ST18M48/" + fileName));
        br = new BufferedReader(new FileReader(fileName));

        readIt();
    }

    private String checkline(String line, int lineNumber, int size, String lineType)
    {
        try
        {
		int temp;
		float floatTemp;
            String[] str = line.split("\t");

		for(int i=0;i<str.length;i++)
			str[i] = str[i].trim();
		//System.out.println(line + "\t" + str.length);

            if(str.length!=size)
                throw new Exception("error");

	    if(lineType.startsWith("S"))
	    {
		//S       00013   00013   1       136     shamu041        1439.88 757.5   18.2    13998110
		//S       00007   00007   1       33      shamu030        382.34  1650.9  0.0     2659620
		if(str[0].length() != 1)
			throw new Exception("error");
		
		temp = Integer.parseInt(str[1]);
		temp = Integer.parseInt(str[2]);
		temp = Integer.parseInt(str[3]);
		temp = Integer.parseInt(str[4]);

		if(str[5].indexOf(" ")>0) throw new Exception("error");

		floatTemp = Float.parseFloat(str[6]);
		floatTemp = Float.parseFloat(str[7]);
		floatTemp = Float.parseFloat(str[8]);
		temp = Integer.parseInt(str[9]);

            } else if(lineType.startsWith("M"))
 	    {
		//M         1       2     1439.690        0.0000   0.2457   20.0    2     20        C.KCKECEKTFHW.S       U
		temp = Integer.parseInt(str[1]);
		temp = Integer.parseInt(str[2]);
		floatTemp = Float.parseFloat(str[3]);
		floatTemp = Float.parseFloat(str[4]);
		floatTemp = Float.parseFloat(str[5]);
		floatTemp = Float.parseFloat(str[6]);
		temp = Integer.parseInt(str[7]);
		temp = Integer.parseInt(str[8]);

/*
		temp = str[9].length()-2;

		for(int i=2;i<temp;i++)
		{
			char ch = str[9].charAt(i);
			if( !Character.isLetter(ch) && ch!='#' && ch!='@' && ch!='*')
				 throw new Exception("error");

		}
*/

		if(str[10].length()!=1 || !Character.isLetter(str[10].charAt(0)) )
			throw new Exception("error");
			
	    }

        }
        catch(Exception e)
        {
            System.out.println("Error=====>> " +  lineType  + "Corrupted : Line # " + lineNumber);

	    StringBuffer sb = new StringBuffer();
	    sb.append("Error=====>> ");
	    sb.append(lineType);
	    sb.append("Corrupted : Line # ");
	    sb.append(lineNumber);
	    sb.append("<BR>"); //for the web display
	    return sb.toString();
        }

	return "";
    }

    private void readIt() throws IOException, Exception
    {
	
        ArrayList errorList = new ArrayList();

	lastLine = br.readLine();
        while( lastLine.startsWith("H\t") || lastLine.startsWith(" ")) {
            header.append(lastLine).append("\n");
            lastLine = br.readLine(); lineNumber++;
        } //read head line

        while( lastLine.equals("") ) { lastLine = br.readLine(); lineNumber++; }//read empty line
       //Check S line

        //System.out.println(lastLine);
	if(!lastLine.startsWith("S"))
            System.out.println("First S line seems to be corrupted");

	//StringBuilder sb = new StringBuilder();
        StringBuffer eachSpec;
        boolean isError = false;

        String result=null;
        StringBuffer fixedData = new StringBuffer();

        while( lastLine!=null && lastLine.startsWith("S") )
        {
            isError = false;
            eachSpec = new StringBuffer();

            String sLine = lastLine;

            result = checkline(lastLine, lineNumber, 10, "S Line ");
            if(!"".equals(result))
                isError = true;

            eachSpec.append( lastLine ).append("\n");

            while ( ( ++lineNumber>0 && (lastLine = br.readLine()) != null))
            {

                //System.out.println(lastLine);
                while(lastLine!=null && lastLine.startsWith("M"))
                {
                    //check M line
                    result = checkline(lastLine, lineNumber, 11, "M Line ");

                    if(!"".equals(result))
                        isError = true;

                    eachSpec.append( lastLine ).append("\n");


                    while ( (++lineNumber>0 && (lastLine = br.readLine()) != null && lastLine.startsWith("L")) || "".equals(lastLine) )
                    {
			if("".equals(lastLine)) continue;

                        result = checkline(lastLine, lineNumber, 2, "L Line ");

                        if(!"".equals(result))
                            isError = true;

                        eachSpec.append( lastLine ).append("\n");
                    }
                }

                if(lastLine==null || lastLine.startsWith("S\t"))
                   break;		
            }

            if(!isError)
                fixedData.append(eachSpec.toString());
            else
                errorList.add(sLine);


		while(lastLine != null && !lastLine.startsWith("S\t"))
		{
			if(lastLine.equals("")) continue;

			if(!"".equals(lastLine))
			{
				System.out.println("Error=====>> Unexpected line : Line  + # " + lineNumber  + "\t" + lastLine);
			}

			lastLine = br.readLine();
			lineNumber++;
		}
        }

        //backup the sqt file
        if(!fixit)
            return;

        if(errorList.size()>0)
        {
            StringBuffer command = new StringBuffer();

            if(backup)
            {
                command.append("mv ");
                command.append(fileName).append(" ");
                command.append(fileName).append(".bak");
            } else
            {
                command.append("rm ");
                command.append(fileName);
            }


            Process p = Runtime.getRuntime().exec(command.toString());
            p.waitFor();

            // Connect print stream to the output stream
            PrintStream out = new PrintStream(new FileOutputStream(fileName));

            out.println(header.toString());
            out.println(fixedData.toString());

            if(backup)
            {
                System.out.print(fileName + " was fixed.  Original file was saved as ");
                System.out.println(fileName + ".bak");
            } else
            {
                System.out.print(fileName + " was fixed.");
            }

            out = new PrintStream( new FileOutputStream("validator.log", true) );

            out.print("Spectra removed in ");
            out.println(fileName);

            for(Iterator itr=errorList.iterator(); itr.hasNext(); )
            {
                out.println(itr.next().toString());
            }

            out.println("\n");
            out.close();
        }

    }

    public static Iterator getAllFiles(String path)
    {
        File f = new File(path);
        return Arrays.asList( f.list() ).iterator();
    }

    public void close() throws IOException
    {
        br.close();
    }

    public String getHeader()
    {
        return header.toString();
    }
}
