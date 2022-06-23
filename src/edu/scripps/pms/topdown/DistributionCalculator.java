

/**
 * @file DistributionCalculator.java
 * This is the source file for edu.scripps.pms.mspid.DistributionCalculator
 * @author Tao Xu
 * @date $Date: 2007/04/12 22:49:56 $
 */
package edu.scripps.pms.topdown;

public class DistributionCalculator {

    public static final double [] COF = { 
        76.18009172947146, -86.50532032941677, 24.01409824083091,
	-1.23173957240155, 0.001208650973866179, -0.000005395239384953
    };
    public static final double STP = 2.5066282746310005;
    public static final double EPS = 0.0000003;
    public static final double FPMIN = 0.000000000000000000000000000001;
    public static final int ITMAX = 100;
    public static final int N = 251;
    //public static final double [] FACTORIAL = new double[N];
    // first dim p, second dim n, third dim k
    private static boolean isLoaded = false;
    public DistributionCalculator() {
    }
    public static double getBinomialSum(double percentage, int n, int k) {
        return calcBinomialSum(percentage, n, k);
    }
    public static double factorial(int n) {
        double result = 1;
        for(int i = 1; i <= n; i++) {
            result *= i;
        }
        return result;
    }
    public static double calcBinomialSum(double p, int n, int k) {
        double prob = 0;
        for(int i = k; i <= n; i++) {
            prob += calcBinomialProbability(p, n, i);
        } 
        //System.out.println("ptrue: " + p + "\tn: " + n + "\tk: " + k + "\tprob: " + prob);
        return prob;
    } 
    public static double calcBinomialProbability(double p, int n, int k) {
        double result = 1;
        double q = 1 - p;
        int nMinusK = n - k;
/*
        result *= FACTORIAL[n];
        result /= FACTORIAL[k];
        result /= FACTORIAL[nMinusK];
*/
        result = nChooseK(n, k);       


        result *= Math.pow(p, k);
        result *= Math.pow(q, nMinusK);
       // System.out.println("\tp: " + p + "\tn: " + n + "\tk: " + k + "\tprob: " + result);
        return result;
    }
    /*
     * hypergeometry calculates the hypergeometrical
     * probability: numTotal elements, numMarked of them marked.
     * numChosen elements are drawn from the pool.
     * what it is the probability that numMarkedChosen elements
     * are marked.
     * numTotal >= numChosen >= numChosen-numMarkedChosen; numTotal >= numMarked >= numMarkedChosen;
     */
    public static double hypergeometry(long numTotal, int numMarked, int numChosen, int numMarkedChosen) {
System.out.println("in hypergeomtry, numTotal: " + numTotal + "\tnumMarked: " + numMarked + "\tnumChosen: " + numChosen + "\tnumMarkedChosen: " + numMarkedChosen);
        double ftemp = 1.0;
        int numUnmarkedChosen = numChosen - numMarkedChosen; // num not marked in the draw
        long numUnmarked = numTotal - numMarked; // num not marked in the pool
        int numMarkedNotChosen = numMarked - numMarkedChosen; // num marked not chosen
        long numNotChosen = numTotal - numChosen; // num in the pool not chosen
        if(numTotal < 0 || numMarked < 0 || numNotChosen < 0 || numMarkedNotChosen < 0) {
            throw new InvalidArgumentException ("WRONG input to HYPERGEOMETRY");
        }
        for(int i = 1; i <= numChosen; i++) {
            if(i <= numMarkedChosen && i <=numUnmarkedChosen) {
                ftemp *= (numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen+i);
                ftemp *= (numMarkedNotChosen+i)/(double)i;
            } else if (i <= numMarkedChosen && i > numUnmarkedChosen) {
                ftemp *= (numMarkedNotChosen + i)/(double)(numNotChosen + i);
            } else if (i > numMarkedChosen && i <= numUnmarkedChosen) {
                ftemp *= (numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen + i);
            } else if (i > numMarkedChosen && i > numUnmarkedChosen) {
                ftemp *= (double)(numNotChosen + i);
            } else {
                System.out.println("WRONG Configuration");
                System.exit (1);
            }
        }
        if(ftemp < 0.0 || ftemp >= Double.MAX_VALUE) {
            double dtemp2 = 1.0;
            ftemp = 0.0;
            for(int i = 1; i <= numChosen; i++) {
                if(i <= numMarkedChosen && i <=numUnmarkedChosen) {
                    dtemp2 = 1.0;
                    dtemp2 *= (numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen+i);
                    dtemp2 *= (numMarkedNotChosen+i)/(double)i;
                    ftemp += Math.log(dtemp2);
                } else if (i <= numMarkedChosen && i > numUnmarkedChosen) {
                    ftemp += Math.log((numMarkedNotChosen + i)/(double)(numNotChosen + i));
                } else if (i > numMarkedChosen && i <= numUnmarkedChosen) {
                    ftemp += Math.log((numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen + i));
                } else if (i > numMarkedChosen && i > numUnmarkedChosen) {
                    ftemp += Math.log(i/(double)(numNotChosen + i));
                }
            }
            ftemp = Math.exp(ftemp);
        } else {
            return ftemp; 
        }
        if(ftemp < 0.0 || ftemp >= Double.MAX_VALUE) {
            System.out.println("Failed Accuracy");
        }
        return ftemp;
    }

    public static double chisquare(double expr[], double theo[], int NBin, 
			   int knstrn) {
	double chsq = 0;
	for(int i = 0; i < NBin; i++) {
	    chsq += (theo[i] - expr[i])*(theo[i] - expr[i])/theo[i];
        }
	int df = NBin - knstrn;
	double prob = gammQ((0.5f*df), (0.5f*chsq));
	return prob;
    }

    /*
     *
     * Returns incomplete gamma function
     * Q(a,x) = 1 - P(a,x)
     */
    static double gammQ(double a, double x) {
	if(x < 0 || x < 0) {
	    throw new InvalidArgumentException("Bad arguments to gammQ");
	}

	if(x < a + 1.f) {
	    return (1.f - gser(a, x));
	} else {
	    return gcf(a,x);
	}
    }

    /*
     * gser returns incomplete gamma function P(a,x)
     * evaluated by its series representation
     * Uses gammaLn
     */
    static double gser(double a, double x) {

	if(x < 0) {
	    throw new InvalidArgumentException("Bad argument in gser");
	}
        if (x == 0) {
	    return 0;
        }
	
	double gln = gammaLn(a);
	double ap = a;
	double sum = 1.f/a;
	double del = sum;

	int i = 0; 
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
	return sum*Math.exp(-x+a*Math.log(x) - gln);
    }

    /*
     *
     * gcf - returns the incomplete gamma function Q(a,x)
     * evaluated by its continued fraction representation
     * Uses gammaLn
     */
    static double gcf(double a, double x) {
	double gln = gammaLn(a);

	double b = x + 1 - a;
	double c = 1/FPMIN;
	double d = 1/b;
	double h = d;
	int i = 0;
	for(i = 1; i <= ITMAX; i++) {
	    double an = -i*(i-a);
	    b += 2;
	    d = an*d + b;
	    if(Math.abs(d) < FPMIN) {
		d = FPMIN;
            }
	    c = b + an/c;
	    if(Math.abs(c) < FPMIN) {
		c = FPMIN;
            }
	    d = 1/d;
	    double del = d*c;
	    h *= del;
	    if(Math.abs(del - 1) < EPS) {
		break;
            }
	}
	if(i == ITMAX) {
	    System.out.println("Need higher ITMAX in gcf");
	    System.exit(1);
	}
	return  Math.exp(-x + a*Math.log(x) - gln)*h;
    }

    /*
     *
     * gammaLn the value ln(GAMMA(xx)) for xx > 0;
     */
    static double gammaLn(double xx) {
	double y = xx;
	double tmp = xx + 5.5;
	tmp = (xx + 0.5)*Math.log(tmp) - tmp;
	double ser = 1.000000000190015;
	for(int j=0; j < COF.length; j++) {
	    y += 1;
	    ser += COF[j]/y;
	}
	return tmp + Math.log(STP*ser/xx);
    }
    public static double nChooseK(int n, int k) {
        double result = 1;
        int nMinusK = n - k;
        int smaller = k;
        int greater = nMinusK;
        if(k > nMinusK) {
            greater = k;
            smaller = nMinusK;
        } 
//System.out.print("n: " + n + "\tk: " + k);
        while(smaller > 0 ) {
            result /= smaller--;
            if(n > greater) {
                result *= n--;
            }
        }
//System.out.print("\t" + result);
        while (n > greater) { 
            result *= n--;
        }
//System.out.println("\t" + result);
        return result;
    }
    public static void main(String [] args) {
        
        System.out.println("Staring...");
        DistributionCalculator dc = new DistributionCalculator();
        System.out.println(dc.getBinomialSum(30, 30, 20));
        System.out.println("Finished");
    }
}




