//package edu.scripps.pms.relex;
package edu.scripps.pms.util.isotope;


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
 * @version 1.0
 */
public class ElementComposition
{
    //private String peptide;
    private static final int ISOTOPE_SIZE = 9;
    private int[] sampleArr = new int[ISOTOPE_SIZE];
    private int[] refArr = new int[ISOTOPE_SIZE];
    private IsotopeTable<String, int[]> isoTable=null;
    private char[] peptide;
    private int start;
    private int length;    
    
    public ElementComposition(String peptideStr, IsotopeTable<String, int[]> isoTable)
    {   
        this(peptideStr.toCharArray(), 0, peptideStr.length(), isoTable);
        //System.out.println("peptide-->>" + peptideStr);
    }
    
    /*
     * length : length of char array to calculate
     */
    public ElementComposition(char[] peptide, int start, int length, IsotopeTable<String, int[]> isoTable)
    {
        this.peptide = peptide;
        this.start = start;
        this.length = length;
        this.isoTable = isoTable;
        
        calculate();
        
        //for(int i=0;i<sampleArr.length;i++)
        //    System.out.println("ele==>>" + sampleArr[i]);
    }

    public int[] getElementSampleArr()
    {
        return sampleArr;
    }

    public int[] getElementRefArr()
    {
        return refArr;
    }

    private void sumSampleElement(int[] arr)
    {
        for(int i=0; i<this.ISOTOPE_SIZE; i++)
        {
//            System.out.println(arr[i]);
            sampleArr[i] += arr[i];
        }
    }

    private void sumRefElement(int[] arr)
    {
        for (int i = 0; i < this.ISOTOPE_SIZE; i++)
        {
            refArr[i] += arr[i];
        }

    }

    private void calculate()
    {
        char c;
        int[] arr;

        for(int i=start; i<length; i++)
        {
            sumSampleElement(isoTable.get("sample" + peptide[i]));
            sumRefElement(isoTable.get("ref" + peptide[i]));
        }
        
        sumSampleElement(isoTable.get("sampleNTERM"));
        sumSampleElement(isoTable.get("sampleCTERM"));

        sumRefElement(isoTable.get("refNTERM"));
        sumRefElement(isoTable.get("refCTERM"));
    }
      
    public void calculateBion()
    {  
        // calculate b ion from y ion
        sampleArr[1] -= 1;  //Hydrogen
        sampleArr[2] -= 1;  //Oxygen      
        
        refArr[1] -= 1;
        refArr[2] -= 1;
    }    
}

            /*
            switch( peptide.charAt(i) )
            {
                case 'A':

                    break;

                case 'R':

                    break;
                case 'N':
                    break;
                case 'D':
                    break;
                case 'C':
                    break;
                case 'Q':
                    break;
                case 'E':
                    break;
                case 'G':
                    break;
                case 'H':
                    break;
                case 'I':
                    break;
                case 'L':
                    break;
                case 'K':
                    break;
                case 'M':
                    break;
                case 'F':
                    break;
                case 'P':
                    break;
                case 'S':
                    break;
                case 'T':
                    break;
                case 'W':
                    break;
                case 'Y':
                    break;
                case 'V':
                    break;
                case 'a':
                    break;
                case 'b':
                    break;
                case 'c':
                    break;

                default :
                    break;
             */
