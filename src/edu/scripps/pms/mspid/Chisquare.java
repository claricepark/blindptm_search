/*
 *
 * NULL hypothesis the experimental data is drawn from the
 * theoretical distribution function
 * experimental data is in expr[]
 * theo - theoretical data
 * knstrn is a parameter, if ==1 it means that the
 *  sum(expr[]) = sum(theo[]);
 */
//import java.io.*;
//import java.util.*;
//import java.text.DecimalFormat;

public class Chisquare {
    static float chisquare(float expr[], float theo[], int NBin, 
			   int knstrn) {
	int j;
	float chsq, df,prob;
	chsq = df = prob=0.f;
	df = (float)(NBin - knstrn);
	chsq = 0.f;
	for(int i = 0; i < NBin; i++)
	    chsq += (float)(theo[i] - expr[i])*
		(float)(theo[i] - expr[i])/
		theo[i];
	prob = gammQ((float)(0.5*df),(float)(0.5*chsq));
	return prob;
    }
    /*
     *
     * Returns incomplete gamma function
     * Q(a,x) = 1 - P(a,x)
     */

    static float gammQ(float a, float x) {
	float gamser, gln;
	gamser = gln = 0.f;
	if(x < 0 || x < 0) {
	    System.out.println("Bad arguments to gammQ");
	    System.exit (1);
	}
	if(x < a + 1.f) {
	    gamser = gser(a,x);
	    return (1.f - gamser);
	}
	else {
	    return gcf(a,x);
	}
    }
    /*
     * gser returns incomplete gamma function P(a,x)
     * evaluated by its series representation
     * Uses gammaLn
     */

    static float gser(float a, float x) {
	double EPS, ap,del,sum,gln;
	int i,ITMAX = 100;
	EPS = 0.0000003;
	ap = del = sum = gln = 0.f;
	gln = gammaLn(a);
	if(x <= 0) {
	    if(x < 0) {
		System.out.println("Bad argument in gser");
		System.exit(1);
	    }
	    return 0.f;
	}
	ap = a;
	sum = 1.f/a;
	del = sum;
	for(i=1; i <= ITMAX; i++) {
	    ap += 1.f;
	    del = del*x/ap;
	    sum += del;
	    if(Math.abs(del) < Math.abs(sum)*EPS)
		break;
	}
	if(i == ITMAX) {
	    System.out.println("a too large, ITMAX to small to converge");
	    System.exit (1);
	}
	return (float)(sum*Math.exp(-x+a*Math.log(x) - gln));
    }
    /*
     *
     * gcf - returns the incomplete gamma function Q(a,x)
     * evaluated by its continued fraction representation
     * Uses gammaLn
     */

    static float gcf(float a, float x) {
	double ap,del,sum,gln,b,c,d,h,an;
	int i,ITMAX = 100;
	double EPS = 0.0000003;
	double FPMIN = 0.000000000000000000000000000001;
	ap = del = sum = gln = 0.f;
	gln = gammaLn(a);
	/*
	System.out.println(Math.exp(gammaLn(4.0f)));
	System.out.println(EPS + "\t" + FPMIN);
	System.out.println(gammaLn(4.0f));
	*/
	b = x + 1. - a;
	c = 1./FPMIN;
	d = 1./b;
	h = d;
	for(i=1; i <= ITMAX; i++) {
	    an = -i*(i-a);
	    b = b+2.;
	    d = an*d + b;
	    if(Math.abs(d) < FPMIN)
		d = FPMIN;
	    c = b + an/c;
	    if(Math.abs(c) < FPMIN)
		c = FPMIN;
	    d = 1./d;
	    del = d*c;
	    h = h*del;
	    if(Math.abs(del - 1.) < EPS)
		break;
	}
	if(i == ITMAX) {
	    System.out.println("Need higher ITMAX in gcf");
	    System.exit(1);
	}
	return  (float)(Math.exp(-x + a*Math.log(x) - gln)*h);
    }
    /*
     *
     * gammaLn the value ln(GAMMA(xx)) for xx > 0;
     */

    static float gammaLn(float xx) {
	int j;
	double ser, stp, tmp,x,y;
	double cof[] = new double [6];
	cof[0] = 76.18009172947146;
	cof[1] = -86.50532032941677;
	cof[2] = 24.01409824083091;
	cof[3] = -1.23173957240155;
	cof[4] = 0.001208650973866179;
	cof[5] = -0.000005395239384953;
	stp = 2.5066282746310005;
	x = xx;
	y = x;
	tmp = x + 5.5d;
	tmp = (x + 0.5d)*Math.log(tmp) - tmp;
	ser = 1.000000000190015;
	for(j=0; j < cof.length; j++) {
	    //System.out.println(cof[j]);
	    y = y + 1.d;
	    ser = ser + cof[j]/y;
	}
	return (float)(tmp + Math.log(stp*ser/x));
    }
}
