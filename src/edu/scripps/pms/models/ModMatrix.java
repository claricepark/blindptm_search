package edu.scripps.pms.models;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class ModMatrix
{
    private final int MATRIX_WIDTH = 13;
    private final int MATRIX_BASIC_HEIGHT = 20;

    private int[][] arr;
    private int[][] rArr;
    private int[][] gArr;
    private int[][] bArr;

    private int width;
    private int height;

    private String residues;

    public ModMatrix(int height, int width, String residues)
    {
        this(height, width);
        this.residues = residues;
    }

    public ModMatrix(int height, int width)
    {
        this.width = width;
        this.height = height;

        arr = new int[height][width];
        rArr = new int[height][width];
        gArr = new int[height][width];
        bArr = new int[height][width];
    }

    public void setArr(int[][] arr)
    {
        this.arr = arr;
    }

    public void setWidth(int width)
    {
        this.width = width;
    }

    public void setHeight(int height)
    {
        this.height = height;
    }

    public int[][] getArr()
    {
        return arr;
    }

    public int getWidth()
    {
        return width;
    }

    public int getHeight()
    {
        return height;
    }

    public void addValue(int i, int j)
    {
        ++arr[i][j];
    }

    public int getValue(int i, int j)
    {
        return arr[i][j];
    }

    public void setValue(int i, int j, int value)
    {
        arr[i][j] = value;
    }

    //maximum non modification value
    public int maxValue()
    {
        int max=0;
        //int temp;

        for(int i=0;i<MATRIX_BASIC_HEIGHT;i++)
        {
            for(int j=0;j<MATRIX_WIDTH;j++)
            {
                if( max<arr[i][j])
                {
                    max = arr[i][j];
                }
            }
        }

        return max;
    }

    //Maximum modification value
    public int maxModValue()
    {
        int max=0;

        for(int i=MATRIX_BASIC_HEIGHT;i<height;i++)
        {
            for(int j=0;j<MATRIX_WIDTH;j++)
            {
                if( max<arr[i][j])
                {
                    max = arr[i][j];
                }
            }
        }

        return max;
    }


    public String getResidues()
    {
        return residues;
    }
}
