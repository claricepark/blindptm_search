package edu.scripps.pms.util.stats;

import java.io.*;
import java.util.StringTokenizer;
import java.text.DecimalFormat;
import java.text.NumberFormat;

public class Fisher
{
    private static DecimalFormat format = new DecimalFormat("0.0000000");

    public Fisher()
    {
    }

    public static void main(String args[]) throws Exception
    {
    /*
        try
        {
            rightTailedFisherFromFile(args[0], args[1], args[2]);
        }
        catch(FileNotFoundException f)
        {
            System.out.println((new StringBuilder()).append("unable  to read ").append(args[0]).toString());
            System.exit(1);
        }
        catch(IOException i)
        {
            System.out.println(i.getMessage());
            System.exit(2);
        }
*/

	BufferedReader br = new BufferedReader(new FileReader(args[0]));
	String eachLine = null;
	Fisher f = new Fisher();
	    //System.out.println( f.calculateFisherTwoTail(4135, 691, 23789, 4373) );
	    //System.out.println( f.getFisherRightTail() );

	while( null != (eachLine=br.readLine()) )
	{
	    String[] arr = eachLine.split("\t");
	    if(arr.length<4)
		continue;

	    for(int i=0;i<arr.length;i++)
	    {
		System.out.print(arr[i] + "\t");
	    }

	    f = new Fisher();
	    //double twoPValue = f.calculateFisherTwoTail(4135, Integer.parseInt(arr[1]), 23789, Integer.parseInt(arr[3]));
	    double twoPValue = f.calculateFisherTwoTail(7641, Integer.parseInt(arr[1]), 42079, Integer.parseInt(arr[3]));
	    //double rightPValue = f.rightTailedFisher(4135, Integer.parseInt(arr[1]), 23789, Integer.parseInt(arr[3]));
	    double rightPValue = f.getFisherRightTail();

	    System.out.print( format.format(twoPValue) );
	    System.out.print("\t");
	    System.out.print( format.format(rightPValue) );
	    System.out.print("\t");

	    if(rightPValue<=0.05)
		System.out.println("*");
	    else
		System.out.println("");
	}
    }

    public static void rightTailedFisherFromFile(String fileName, String theOverallTotal, String theOverallChanged)
        throws IOException
    {
        FileReader fr = new FileReader(fileName);
        BufferedReader bIn = new BufferedReader(fr);
        Fisher f = new Fisher();
        int overallTotal = Integer.valueOf(theOverallTotal.trim()).intValue();
        int overallChanged = Integer.valueOf(theOverallChanged.trim()).intValue();
        String line;
        while((line = bIn.readLine()) != null) 
        {
            StringTokenizer tokenizer = new StringTokenizer(line);
            int n4 = Integer.parseInt(tokenizer.nextToken());
            int n2 = Integer.parseInt(tokenizer.nextToken());
            int n = Integer.parseInt(tokenizer.nextToken());
            f.calculateFisherTwoTail(overallChanged, n2, overallTotal, n4);
            System.out.println((new StringBuilder()).append(n4).append("\t").append(n2).append('\t').append(f.getFisherRightTail()).append('\t').append(n).toString());
        }
    }

    public double rightTailedFisher(int theOverallChanged, int theChangedInNode, int theOverallTotal, int theTotalInNode)
    {
        calculateFisherTwoTail(theOverallChanged, theChangedInNode, theOverallTotal, theTotalInNode);
        return Math.log(getFisherRightTail()) / Math.log(10D);
    }

    public double calculateFisherTwoTail(int totalChanged, int changedInNode, int total, int inNode)
    {
        return calculateFisherFromMatrix(changedInNode, totalChanged - changedInNode, inNode - changedInNode, total - inNode - (totalChanged - changedInNode));
    }

    private double calculateFisherFromMatrix(int n11, int n12, int n21, int n22)
    {
        n11_ = n11;
        n12_ = n12;
        n21_ = n21;
        n22_ = n22;
        if(n11_ < 0)
            n11_ *= -1;
        if(n12_ < 0)
            n12_ *= -1;
        if(n21_ < 0)
            n21_ *= -1;
        if(n22_ < 0)
            n22_ *= -1;
        int n1_ = n11_ + n12_;
        int n_1 = n11_ + n21_;
        int n = n11_ + n12_ + n21_ + n22_;
        exact(n11_, n1_, n_1, n);
        left = sless;
        right = slarg;
        twotail = sleft + sright;
        if(twotail > 1.0D)
            twotail = 1.0D;
        return twotail;
    }

    public double getFisherLeftTail()
    {
        return left;
    }

    public double calculateFisherLeftTail(int totalChanged, int changedInNode, int total, int inNode)
    {
        calculateFisherFromMatrix(changedInNode, totalChanged - changedInNode, inNode - changedInNode, total - inNode - (totalChanged - changedInNode));
        return getFisherLeftTail();
    }

    public double getFisherRightTail()
    {
        return right;
    }

    public double calculateFisherRightTail(int totalChanged, int changedInNode, int total, int inNode)
    {
        calculateFisherFromMatrix(changedInNode, totalChanged - changedInNode, inNode - changedInNode, total - inNode - (totalChanged - changedInNode));
        return getFisherRightTail();
    }

    public double calculateFisherTwoTail()
    {
        return twotail;
    }

    private static double lngamm(int z)
    {
        double x = 0.0D;
        x += 1.6594701874084621E-007D / (double)(z + 7);
        x += 9.9349371139307475E-006D / (double)(z + 6);
        x -= 0.1385710331296526D / (double)(z + 5);
        x += 12.50734324009056D / (double)(z + 4);
        x -= 176.61502914983859D / (double)(z + 3);
        x += 771.32342877576741D / (double)(z + 2);
        x -= 1259.1392167222889D / (double)(z + 1);
        x += 676.52036812188351D / (double)z;
        x += 0.99999999999951827D;
        return (Math.log(x) - 5.5810614667953278D - (double)z) + ((double)z - 0.5D) * Math.log((double)z + 6.5D);
    }

    private static double lnfact(int n)
    {
        if(n <= 1)
            return 0.0D;
        else
            return lngamm(n + 1);
    }

    private static double lnbico(int n, int k)
    {
        return lnfact(n) - lnfact(k) - lnfact(n - k);
    }

    private static double hyper_323(int n11, int n1_, int n_1, int n)
    {
        return Math.exp((lnbico(n1_, n11) + lnbico(n - n1_, n_1 - n11)) - lnbico(n, n_1));
    }

    private double hyper(int n11)
    {
        return hyper0(n11, 0, 0, 0);
    }

    private double hyper0(int n11i, int n1_i, int n_1i, int ni)
    {
        if((n1_i == 0) & (n_1i == 0) & (ni == 0))
        {
            if(n11i % 10 != 0)
            {
                if(n11i == sn11 + 1)
                {
                    sprob = sprob * (((double)sn1_ - (double)sn11) / (double)n11i) * (((double)sn_1 - (double)sn11) / (((double)n11i + (double)sn) - (double)sn1_ - (double)sn_1));
                    sn11 = n11i;
                    return sprob;
                }
                if(n11i == sn11 - 1)
                {
                    sprob = sprob * ((double)sn11 / ((double)sn1_ - (double)n11i)) * ((((double)sn11 + (double)sn) - (double)sn1_ - (double)sn_1) / ((double)sn_1 - (double)n11i));
                    sn11 = n11i;
                    return sprob;
                }
            }
            sn11 = n11i;
        } else
        {
            sn11 = n11i;
            sn1_ = n1_i;
            sn_1 = n_1i;
            sn = ni;
        }
        sprob = hyper_323(sn11, sn1_, sn_1, sn);
        return sprob;
    }

    private double exact(int n11, int n1_, int n_1, int n)
    {
        int max = n1_;
        if(n_1 < max)
            max = n_1;
        int min = (n1_ + n_1) - n;
        if(min < 0)
            min = 0;
        if(min == max)
        {
            sless = 1.0D;
            sright = 1.0D;
            sleft = 1.0D;
            slarg = 1.0D;
            return 1.0D;
        }
        double prob = hyper0(n11, n1_, n_1, n);
        sleft = 0.0D;
        double p = hyper(min);
        int i;
        for(i = min + 1; p < 0.99999998999999995D * prob; i++)
        {
            sleft += p;
            p = hyper(i);
        }

        i--;
        if(p < 1.0000000099999999D * prob)
            sleft += p;
        else
            i--;
        sright = 0.0D;
        p = hyper(max);
        int j;
        for(j = max - 1; p < 0.99999998999999995D * prob; j--)
        {
            sright += p;
            p = hyper(j);
        }

        j++;
        if(p < 1.0000000099999999D * prob)
            sright += p;
        else
            j++;
        if(Math.abs(i - n11) < Math.abs(j - n11))
        {
            sless = sleft;
            slarg = (1.0D - sleft) + prob;
        } else
        {
            sless = (1.0D - sright) + prob;
            slarg = sright;
        }
        return prob;
    }

    private double left;
    private double right;
    private double twotail;
    private double sleft;
    private double sright;
    private double sless;
    private double slarg;
    private int sn11;
    private int sn1_;
    private int sn_1;
    private int sn;
    private double sprob;
    private int n11_;
    private int n12_;
    private int n21_;
    private int n22_;
}
