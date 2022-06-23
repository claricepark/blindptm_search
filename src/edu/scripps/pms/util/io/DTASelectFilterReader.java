package edu.scripps.pms.util.io;

import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;

import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.ProteinGroup;

/**
 * @author  Robin Park
 * @version $Id: DTASelectFilterReader.java,v 1.29 2012/03/11 06:53:04 taoxu Exp $
 */
public class DTASelectFilterReader
{
    private final long READ_FROM_THE_END = 500; //position from the end

    private String dbFileName;
    private String dbFilePathAndName;
    //    private InputStreamReader reader;
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private String criteria;
    private int unfilteredProteinNum;
    private int redundantProteinNum;
    private int nonRedundantProteinNum;
    private int unfilteredPeptideNum;
    private int redundantPeptideNum;
    private int nonRedundantPeptideNum;

    private boolean version2;  //DTASelect version

    //read total peptide number
    public int getTotalPeptideNumber() throws IOException
    {
	int totalPeptideCount=0;

	DTASelectFilterReader dtaReader = new DTASelectFilterReader(fileName);

	for (Iterator<Protein> itr1 = dtaReader.getProteins(); itr1.hasNext(); ) {
	    Protein protein = itr1.next();

	    totalPeptideCount += protein.getPeptideSize();
	}

	return totalPeptideCount;

    }

    public static void main(String args[]) throws IOException
    {
	DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );

	Iterator<Protein> pitr = reader.getProteins();

	for (Iterator<Protein> itr = pitr; itr.hasNext(); )
	{
	    Protein protein = itr.next();
	    System.out.println(protein.getLocus());

	}

    }

    public DTASelectFilterReader(String fileName) throws IOException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        init();
    }

    /*
    public Iterator <Protein> getProteins(final InputStream is) throws IOException {
        return getProteins(new InputStreamReader(is));
    }
    */
    private void readSummary() throws IOException
    {
        RandomAccessFile file = new RandomAccessFile(fileName, "r");

        file.seek(file.length()-500);

        String eachLine;
        eachLine=file.readLine();

        while( (eachLine=file.readLine()) != null && !eachLine.startsWith("Unfiltered") );

        String[] arr = eachLine.split("\t");
        this.unfilteredProteinNum = Integer.parseInt(arr[1]);
        this.unfilteredPeptideNum = Integer.parseInt(arr[2]);

        eachLine=file.readLine(); //Redundant line of DTASelect-filter.txt
        arr = eachLine.split("\t");
        this.redundantProteinNum = Integer.parseInt(arr[1]);
        this.redundantPeptideNum = Integer.parseInt(arr[2]);

        eachLine=file.readLine(); //NonRedundant line of DTASelect-filter.txt
        arr = eachLine.split("\t");
        this.nonRedundantProteinNum = Integer.parseInt(arr[1]);
        this.nonRedundantPeptideNum = Integer.parseInt(arr[2]);

        file.close();
    }

    private void init() throws IOException
    {
//        readSummary();

//System.out.println("filename == >>"  + fileName);

        br = new BufferedReader(new FileReader(fileName));

        readHeader();

        //Move line to parameters.  We assume parameters start after carriage return
        while ((lastLine = br.readLine()) != null) {
            if(lastLine.startsWith("ProLuCID") || lastLine.startsWith("SEQUEST") || lastLine.startsWith("BlindPTM")) {
                break;
            }
        }

        StringBuffer sb = new StringBuffer();
        sb.append(lastLine=br.readLine());  //Read this line, which can be either empty or not

        while (!(lastLine = br.readLine()).equals(""))
        {
            sb.append(lastLine);
            sb.append("\n");
        }

        criteria = sb.toString();

        //Move line to parse protein
        while (!(lastLine = br.readLine()).startsWith("Unique"));
        //add setFeatureIndices here
        Peptide.setFeatureIndices(lastLine);
        lastLine = br.readLine();
    }

    private void readHeader() throws IOException
    {
        lastLine = br.readLine();

        if(lastLine.startsWith("DTASelect v2.0"))
            version2 = true;
        else
            version2 = false;

//        System.out.println(lastLine.startsWith("DTASelect v2.0"));
//        System.out.println(lastLine);

        //whie (!(lastLine = br.readLine()).endsWith("fasta")); .startsWith("Unique"));
        for(int i=0; i<2; i++)
            lastLine = br.readLine();

        //Remove directory name
        this.dbFilePathAndName = lastLine;
        this.dbFileName = lastLine.substring(lastLine.lastIndexOf("/")+1);
    }
    public ArrayList<Protein> getProteinList() throws IOException {
        ArrayList<Protein> list = new ArrayList();
        Iterator <Protein> it = getProteins(); 
        while(it.hasNext()) {
            list.add(it.next());
        } 
        return list;
    }
    
    public ArrayList<ProteinGroup> getProteinGroupList() throws IOException {
        ArrayList<ProteinGroup> list = new ArrayList();
        Iterator <Protein> it = getProteins(); 
        while(it.hasNext()) {
            Protein p = it.next();
            if(p.getNumPeptides() == 0) {
                ProteinGroup pg = new ProteinGroup();
                pg.addProtein(p);
                list.add(pg);
                while(it.hasNext()) {
                    p = it.next();
                    pg.addProtein(p);
                    if(p.getNumPeptides() > 0) {
                        break;
                    }
                }
              
            } else {
                ProteinGroup pg = new ProteinGroup();
                pg.addProtein(p);
                list.add(pg);
            }      
        } 

        return list;
    }
    public Iterator <Protein> getProteins() throws IOException {
        return new Iterator<Protein>() {
            private Protein protein;
            private Peptide peptide;

            public boolean hasNext() {
                return lastLine != null && !lastLine.startsWith("\tProteins\t");
            }

            public Protein next() {

                try {
                    protein = getProtein(lastLine.split("\t"));
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private Protein getProtein(String[] strArr) throws IOException {

//                String[] strArr = lastLine.split("\t");
//		String[] peptideLine;
                //protein = new Protein(strArr);
                protein = new Protein();
                try {
                    protein.setElement(strArr);
                } catch (ArrayIndexOutOfBoundsException e ) {
                    System.out.println(e.getMessage());
                }
                /*
                The format of the DTASelect-filter.txt file peptide line
                 is the following:

                 - a star if the peptide is unique to that protein (optional)
                 - then an integer greater than 1 (optional)
                 - then one of the characters 'M', 'Y', or 'N' (optional)
                 - then tab (mandatory)
                 - then the rest of the fields...

                 Anything that comes until the first tab is optional. You can
                 have a line that starts with "\t", or "*\t", or "2\t", or
                 "*2\t", or "*2M\t", or "Y\t", or "2N\t", or "*Y\t", etc...

		Peptide line starts like (*)(int)(M||Y||N)\t
                 */

                /*
                 *  Some proteins does not have peptide lines, because those proteins
                 *  are assumed to have identical peptides as following protein has.
                 *
                 **/


                while ( ((lastLine = br.readLine()) != null && !lastLine.equals(""))
                    && !lastLine.startsWith("\tProteins\t"))
                {
                    strArr = lastLine.split("\t");

                    // If Spectrum Count position does not have a decimal point,
                    // it is not a peptide line
                    if(strArr[2].indexOf(".")<=0)
                            break;

                    peptide = new Peptide(strArr, version2);
                    protein.addPeptide(peptide);
                }

                return protein;
            }
        };
    }

    public void close() throws IOException
    {
        br.close();
    }

    public void setDbFileName(String dbFileName)
    {
        this.dbFileName = dbFileName;
    }

    public void setUnfilteredProteinNum(int unfilteredProteinNum)
    {
        this.unfilteredProteinNum = unfilteredProteinNum;
    }

    public void setRedundantProteinNum(int redundantProteinNum)
    {
        this.redundantProteinNum = redundantProteinNum;
    }

    public void setNonRedundantProteinNum(int nonRedundantProteinNum)
    {
        this.nonRedundantProteinNum = nonRedundantProteinNum;
    }

    public void setUnfilteredPeptideNum(int unfilteredPeptideNum)
    {
        this.unfilteredPeptideNum = unfilteredPeptideNum;
    }

    public void setRedundantPeptideNum(int redundantPeptideNum)
    {
        this.redundantPeptideNum = redundantPeptideNum;
    }

    public void setNonRedundantPeptideNum(int nonRedundantPeptideNum)
    {
        this.nonRedundantPeptideNum = nonRedundantPeptideNum;
    }

    public String getDbFileName()
    {
        return dbFileName;
    }
    public String getDbFilePathAndName()
    {
        return dbFilePathAndName;
    }

    public String getCriteria()
    {
        return criteria;
    }

    public int getUnfilteredProteinNum()
    {
        return unfilteredProteinNum;
    }

    public int getRedundantProteinNum()
    {
        return redundantProteinNum;
    }

    public int getNonRedundantProteinNum()
    {
        return nonRedundantProteinNum;
    }

    public int getUnfilteredPeptideNum()
    {
        return unfilteredPeptideNum;
    }

    public int getRedundantPeptideNum()
    {
        return redundantPeptideNum;
    }

    public int getNonRedundantPeptideNum()
    {
        return nonRedundantPeptideNum;
    }

    public boolean isVersion2() {
        return version2;
    }
}
