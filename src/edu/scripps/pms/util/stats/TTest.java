package edu.scripps.pms.util.stats;

/**
 *
 * @author  Robin Park
 * @version $Id: TTest.java,v 1.3 2007/03/09 18:21:09 taoxu Exp $
*/

public class TTest 
{
    public static void main(String args[]) throws Exception
    {
	TTest dev = new TTest();

	System.out.println( dev.T_p(1.96, 7) );
        double z = 1;
        System.out.println(z + "\t" + Norm_p(z));
        z = 2;
        System.out.println(z + "\t" + Norm_p(z));
        z = 3;
        System.out.println(z + "\t" + Norm_p(z));
        z = 4;
        System.out.println(z + "\t" + Norm_p(z));
        z = 1.96;
        System.out.println(z + "\t" + Norm_p(z));
        z = 2.36;
        System.out.println(z + "\t" + Norm_p(z));
    }

    public static double T_p(double t, double df) {
	double abst = Math.abs(t);
	double tsq = t*t;
	double p;
	    //var abst = abs(t), tsq = t*t, p;
	if(df == 1)
	    p = 1 - 2*Math.atan(abst)/Math.PI;
	else if(df == 2)
	    p = 1 - abst/Math.sqrt(tsq + 2);
	else if(df == 3) 
	    p = 1 - 2*(Math.atan(abst/Math.sqrt(3)) + abst*Math.sqrt(3)/(tsq + 3))/Math.PI;
	else if(df == 4) 
	    p = 1 - abst*(1 + 2/(tsq + 4))/Math.sqrt(tsq + 4);
	else 
	{
	    // finds the z equivalent of t and df st they yield same probs.
	    double z = T_z(abst, df);
	    if (df>4) 
		p = Norm_p(z);
	    else 
		p = Norm_p(z); // small non-integer df
	}

	return p;
    }

    public static double Norm_p(double z) 
    {
	// Returns the two-tailed standard normal prob. level given a z.
	double absz = Math.abs(z);
	double a1 = 0.0000053830, a2 = 0.0000488906, a3 = 0.0000380036,
	    a4 = 0.0032776263, a5 = 0.0211410061, a6 = 0.0498673470;
	double p = (((((a1*absz+a2)*absz+a3)*absz+a4)*absz+a5)*absz+a6)*absz+1;
	p = Math.pow(p, -16);
	return p;
    }

    public static double T_z(double t, double df) {
	// Converts a t value to an approximate z value w.r.t the given df
	// s.t. std.norm.(z) = t(z, df) at the two-tail probability level.

	double A9 = df - 0.5;
	double B9 = 48*A9*A9;
        double T9 = t*t/df, Z8, P7, B7, z;

	if (T9 >= 0.04) 
	    Z8 = A9*Math.log(1+T9);
	else  
	    Z8 = A9*(((1 - T9*0.75)*T9/3 - 0.5)*T9 + 1)*T9;
	P7 = ((0.4*Z8 + 3.3)*Z8 + 24)*Z8 + 85.5;
	B7 = 0.8*Math.pow(Z8, 2) + 100 + B9;
	z = (1 + (-P7/B7 + Z8 + 3)/B9)*Math.sqrt(Z8);

	return z;
    }

}
