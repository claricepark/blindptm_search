package edu.scripps.pms.util.seq;

import edu.scripps.pms.helper.PMSConstants;
import java.io.*;
import java.util.*;
import edu.scripps.pms.models.DbSearchParameters;
import edu.scripps.pms.util.dtaselect.ModResidue;
import gnu.trove.TIntDoubleHashMap;

/**
 * @author  Robin Park
 * @version $Id: ParamReader.java,v 1.15 2007/03/05 23:09:22 rpark Exp $
 */
public class ParamReader
{
    private String path;
    private String fileName;
    private DbSearchParameters param;
    private boolean isModification;
    private Hashtable<String, String> ht = new Hashtable<String, String> ();
    private List residueList = new Vector();

    private final char[] MOD_SYMBOL = {'*', '#', '@', };
    private TIntDoubleHashMap symbolMap = new TIntDoubleHashMap(10);

    private double[] massShiftArr = new double[3];

    public static void main(String args[]) throws Exception
    {
        ParamReader p = new ParamReader(args[0], args[1], Integer.parseInt(args[2]));
        Hashtable<String, String> h = p.getHashtable();

        System.out.println(h.get("diff_search_options"));
        System.out.println( p.isModification() );
        List l = p.getResidueList();

        for(int i=0;i<l.size();i++)
        {
            ModResidue residue = (ModResidue)l.get(i);

            System.out.println(residue.getMassShift() + "\t" + residue.getResidue());
        }

        double[] d = p.getMassShiftArr();
        for(int i=0;i<d.length;i++)
        {
            System.out.println(d[i]);
        }
    }

    // searchType = 1 for sequest
    public ParamReader(String path, String fileName, int searchType) throws IOException
    {
        this.path = path;
        this.fileName = fileName;

        switch(searchType)
        {
            case PMSConstants.SEQUEST:
                sequestInit();
                break;

            case PMSConstants.PEP_PROBE:
                pepProbeInit();
                break;

            case PMSConstants.GUTENTAG:
                gutentagInit();
                break;
        }
    }

    public void sequestInit() throws IOException
    {
        path = path + File.separator + fileName;
        BufferedReader br = null;
        try
        {
            br = new BufferedReader(new FileReader(path));
            String eachLine;

            //Hashtable<String> ht = new Hashtable<String>();

            String[] strArr;
            StringBuffer sb = new StringBuffer();
            //br.skip(2000);
            while ( (eachLine = br.readLine()) != null)
            {
                sb.append(eachLine);
                sb.append("\n");

                if (eachLine.startsWith("#") ||
                    (strArr = eachLine.split("=")).length < 2)
                    continue;

                ht.put(strArr[0].trim(), strArr[1]);
            }

            param = new DbSearchParameters();
            param.setDb( trimValue(ht.get("database_name")) );
            param.setParametersFile(path);
            param.setParameters(sb.toString());

	    String pepTolerance = ht.get("peptide_mass_tolerance");
	    if(null == pepTolerance)
		pepTolerance = ht.get("ppm_peptide_mass_tolerance");
	   
            param.setPeptideMassTolerance( trimValueAsFloat(pepTolerance) );
            param.setFragmentIonTolerance( trimValueAsFloat(ht.get("fragment_ion_tolerance")) );
            param.setMassTypeParent( trimValueAsInt(ht.get("mass_type_parent")) );
            param.setMassTypeFragment( trimValueAsInt(ht.get("mass_type_fragment")) );
            param.setEnzymeNumber( trimValueAsInt(ht.get("enzyme_number")) );

            if(null != ht.get("diff_search_options") )
            {
                String[] modArr = ht.get("diff_search_options").toString().trim().split(" ");

                ModResidue residue;
                //Read each set of residue
                //e.g. diff_search_options = 80.0 ST -18.0 ST 0.0 X
                //*, #, @

                for(int i=0;i<5;i+=2)
                {
                    double massShift = Double.parseDouble(modArr[i]);

                    if( massShift != 0 )
                    {
                        symbolMap.put( MOD_SYMBOL[i/2], massShift);

                        this.isModification = true;

                        for(int j=0;j<modArr[i+1].length();j++)
                        {
                            residue = new ModResidue(modArr[i+1].charAt(j), massShift);
                            massShiftArr[i/2] = massShift;
                            residueList.add(residue);

                        }
                    }
                }
            }

        }
        catch(IOException e)
        {
            System.out.println("Error reading param file " + e);
            throw new IOException(e.toString());
        }
        finally
        {
            br.close();
        }
    }

    public double[] getMassShiftArr()
    {
        return massShiftArr;
    }

    public void pepProbeInit()
    {
    }

    public void gutentagInit()
    {

    }

    public DbSearchParameters getDbSearchParameters()
    {
        return param;
    }

    public int trimValueAsInt(String str)
    {
        return Integer.valueOf( trimValue(str) );
    }

    public float trimValueAsFloat(String str)
    {
        return Float.valueOf( trimValue(str) );
    }

    public String trimValue(String str)
    {
        int index = str.indexOf(';');
        if (index > 0)
            str = str.substring(0, index);

        return str.trim();
    }

    public boolean isModification() {
        return isModification;
    }

    public Hashtable<String, String> getHashtable()
    {
        return ht;
    }

    public List getResidueList() {
        return residueList;
    }

    public TIntDoubleHashMap getSymbolMap()
    {
        return symbolMap;
    }
}
