package edu.scripps.pms.tools;

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

public class IsotopeDistribution
{
    private static final double abund[][] = new double[9][4];
    private static final int DIST_SIZE = 200;
    private static int npeak[] = new int[9];
    private double startMass;
    private double endMass;
    private double avgMass;
    private double[] masslist;
    private double[] relabun = new double[100];

    private int[] element;
    private static double c13=1.003354;

    private final double ISODIST_THRESHOLD=0.10; // remove masses below 10%
    
//    private static double enrichment = 0.98;
    private static double enrichment;
    //= 0.98;
    public void init()//static
    {
//	enrichment = 0.94;
        // Carbon
        abund[0][0] = 100.0;
        abund[0][1] = 1.0958793;
        npeak[0] = 2;

        // Hydrogen
        abund[1][0] = 100.0;
        abund[1][1] = 0.014502102;
        npeak[1] = 2;

        // Oxygen
        abund[2][0] = 100.0;
        abund[2][1] = 0.03799194;
        abund[2][2] = 0.20499609;
        npeak[2] = 3;

        // Nitrogen
        abund[3][0] = 100.0;
        abund[3][1] = 0.368351851;
        npeak[3] = 2;

        // Sulfur
        abund[4][0] = 100.0;
        abund[4][1] = 0.784;
        abund[4][2] = 4.442;
        abund[4][3] = 0.014;
        npeak[4] = 4;

        // Phosphorous
        abund[5][0] = 100.0;
        npeak[5] = 4;

        // 15N Nitrogens
        abund[6][0] = 100 - enrichment * 100;
        abund[6][1] = (enrichment * 100) + (abund[5][0] * abund[3][1] / 100);
        npeak[6] = 2;

        // 2H Deuterium
        abund[7][0] = 100 - 98.27;
        abund[7][1] = 98.27 + (abund[6][0] * abund[1][1] / 100);
        npeak[7] = 2;

        // 13C Carbons
        abund[8][0] = 100 - enrichment * 100;
        abund[8][1] = (enrichment * 100) + (abund[8][0] * abund[0][1] / 100);
        npeak[8] = 2;
    }


    public IsotopeDistribution(int[] element, boolean isLight, double enrichment)
    {
        this.element = element;
	this.enrichment = enrichment;
	init();

	if(isLight)
	    calculateLight();
	else
	    calculateHeavy();
    }

    public IsotopeDistribution(int[] element, boolean isLight)
    {
        this.element = element;
	
	if(isLight)
	    calculateLight();
	else
	    calculateHeavy();
    }

    public static void main(String args[])
    {
        //int[] element = {77, 129, 27, 23, 1, 0, 0, 0, 0, };
        //int[] element = {77, 129, 27, 0, 1, 0, 23, 0, 0, };
        //int[] element = {4, 8, 2, 1, 0, 5, 0, 0, 0, };
        //int[] element = {4, 9, 5, 1, 0, 1, 0, 0, 0, };
        //int[] element = {57, 99, 22, 17, 0, 0, 0, 0, 0, };
        int[] element = {57, 99, 22, 0, 0, 0, 17, 0, 0, };
        IsotopeDistribution iso = new IsotopeDistribution(element, false);
	//iso.setEnrichment(0.94);


        //System.out.println(iso.getEachAbund(0, 0));

        //long start = System.currentTimeMillis();

        System.out.println( iso.getStartMass() + " " + iso.getEndMass());
    }

    private void calculateLight()
    {
        int i0, i1, j, k, p, q, l, ii, cnt, i;
        double[] CPATT = new double[200];
        double[] D = new double[DIST_SIZE];
        double mass = 0.0;

        double max, prec;
        double thres = 0.05, sum, maxmass, dif;
        double[] masses = new double[6];
        double[] fracabun = new double[100];
        masslist = new double[100];

        p = 0;
        q = 0;
        prec = 0.00001;
        CPATT[0] = 1.0;
        for (j = 0; j <= 8; ++j) {

        int count=1;
            if (element[j] > 0) {

                for (ii = 0; ii < element[j]; ++ii) {

                    for (i1 = 0; i1 <DIST_SIZE; ++i1) D[i1] = 0.0; /* Erase D[] */
                    //D = new double[200]; //this slows down

                    /* Calculate Isotope Distribution */
                    for (k = p; k <= q; ++k) {
                        count++;

                        for (l = 0; l < npeak[j]; ++l) {
                            cnt = k + l;
                            D[cnt] = D[cnt] + CPATT[k] * abund[j][l];
                        }
                    }

                    /* Normalize Intensities */
                    q = q + npeak[j] - 1;

                    max = 0.0;
                    for (k = p; k <= q; ++k) {
                        if (D[k] > max) {
                            max = D[k];
                        }
                    }

                    for (k = p; k <= q; ++k) {
                        D[k] /= max;
                    }

                    /* Eliminate small peaks to the left */
                    for (k = p; k <= q; ++k) {
                        if (D[k] > prec) {
                            p = k;
                            //k = q;
                            break;
                        }
                    }
                    /* Eliminate small peaks to the right */
                    for (k = q; D[k] < prec; --k) {
                        q = k;
                        k = k - 1;
                    }

                    /* Create new isotope pattern */
                    for (i0 = 0; i0 < 200; ++i0) CPATT[i0] = 0.0;
                        /* Clear CPATT[] */

                    for (k = p; k <= q; ++k) {
                        CPATT[k] = D[k];
                    }
                }

            }

        }

        /* Calculate mono-isotopic mass */
        mass = (element[0] * 12) + (element[1] * 1.007825) +
            (element[2] * 15.9949146)
            + (element[3] * 14.003074) + (element[4] * 31.9720718) +
            (element[5] * 30.9737620) +
            (element[6] * 14.003074 + (element[7] * 1.007825) +
             (element[8] * 12.000000000));

        masses[0] = mass;

        /* Calculate average and fractional masses*/
        sum = 0;
        maxmass = 0;
        dif = 0;
        max = 0;
        for (k = p; k <= q; k++) {
            sum = sum + CPATT[k];
            if (CPATT[k] > max) {
                max = CPATT[k];
                dif = k;
            }
        }

        i = 0;
        avgMass = 0;

        for (k = p; k <= q; k++) {
            relabun[i] = CPATT[k]; //calculates relative abundances
            fracabun[i] = CPATT[k] / sum; //calculates fractional abundances
            masslist[i] = mass + c13*k; //masses of isotopes
            //masslist[i] = mass + k; //masses of isotopes
            avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
            i = i + 1;
        }
        //masses[1] = avgMass; //avg. mass
        //masses[2] = mass + dif; //calculate mass at maximum intensity
        //masses[3] = dif - p + 1; //calculate which ion has the max intensity

        //calculates beginning mass for integration
        ii = 0;

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > thres) {
                //beginmass = masslist[ii] - 0.5;
                this.startMass = masslist[ii];
                break;
            }
        }
//        masses[4] = beginmass;
        //this.startMass = beginmass;


        for (ii = 0; ii <= q; ii++) {
            
            if (relabun[ii] > thres) {                
                this.endMass = masslist[ii];
            }
        }

/*
        System.out.println("start mass==>> " + startMass + " " + endMass);
        
	for(int u=0;u<masslist.length;u++)
	{
	    System.out.println(u + "mass==>>\t" + masslist[u] + "\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}

	for(int u=0;u<relabun.length;u++)
	{
	    System.out.println(u + "abun==>>\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}
//        masses[5] = endmass;
        this.endMass = endmass;
//        System.out.println(endmass);

        /*
System.out.println("masslist");
        for (int ii0 = 0; ii0 < masslist.length; ii0++) {
            System.out.println(masslist[ii0]);

        }


        for (int ii0 = 0; ii0 < masses.length; ii0++) {
            System.out.println("mass " + ii0 + "\t" + masses[ii0]);

        }
*/
    }

    private void calculateHeavy()
    {
//        this.enrichment = conf.getEnrichment();
        int i0, i1, j, k, p, q, l, ii, cnt, i;
        double[] CPATT = new double[200];
        double[] D = new double[DIST_SIZE];
        double mass = 0.0;

        double max, prec;
        double thres = 0.1, sum, maxmass, dif;
        double[] masses = new double[6];
        double[] fracabun = new double[100];
        masslist = new double[100];

        p = 0;
        q = 0;
        prec = 0.00001;
        CPATT[0] = 1.0;
        for (j = 0; j <= 8; ++j) {

        int count=1;
            if (element[j] > 0) {

                for (ii = 0; ii < element[j]; ++ii) {

                    for (i1 = 0; i1 <DIST_SIZE; ++i1) D[i1] = 0.0; /* Erase D[] */
                    //D = new double[200]; //this slows down

                    /* Calculate Isotope Distribution */
                    for (k = p; k <= q; ++k) {
                        count++;

                        for (l = 0; l < npeak[j]; ++l) {
                            cnt = k + l;
                            D[cnt] = D[cnt] + CPATT[k] * abund[j][l];
                        }
                    }

                    /* Normalize Intensities */
                    q = q + npeak[j] - 1;

                    max = 0.0;
                    for (k = p; k <= q; ++k) {
                        if (D[k] > max) {
                            max = D[k];
                        }
                    }

                    for (k = p; k <= q; ++k) {
                        D[k] /= max;
                    }

                    /* Eliminate small peaks to the left */
                    for (k = p; k <= q; ++k) {
                        if (D[k] > prec) {
                            p = k;
                            //k = q;
                            break;
                        }
                    }
                    /* Eliminate small peaks to the right */
                    for (k = q; D[k] < prec; --k) {
                        q = k;
                        k = k - 1;
                    }

                    /* Create new isotope pattern */
                    for (i0 = 0; i0 < 200; ++i0) CPATT[i0] = 0.0;
                        /* Clear CPATT[] */

                    for (k = p; k <= q; ++k) {
                        CPATT[k] = D[k];
                    }
                }

            }

        }


/*
System.out.println("start=====================");
        for(int cpi=0; cpi<D.length; cpi++)
        {
        //	System.out.println(D[cpi]);
        }
    System.out.println("end=====================");
 */
        /* Calculate mono-isotopic mass */
        mass = (element[0] * 12) + (element[1] * 1.007825) +
            (element[2] * 15.9949146)
            + (element[3] * 14.003074) + (element[4] * 31.9720718) +
            (element[5] * 30.9737620) +
            //(element[6] * 14.003074 + (element[7] * 1.007825) +
            (element[6] * 15.00010897 + (element[7] * 1.007825) +
             (element[8] * 12.000000000));

        masses[0] = mass;

        /* Calculate average and fractional masses*/
        sum = 0;
        maxmass = 0;
        dif = 0;
        max = 0;
        for (k = p; k <= q; k++) {
            sum = sum + CPATT[k];
            if (CPATT[k] > max) {
                max = CPATT[k];
                dif = k;
            }
        }

        i = 0;
        avgMass = 0;

//	double N15=0.997034;
//	mass += p*N15;
	//element[6] = # of nitrogens
        for (k = p; k <= q; k++) {
//	int tempQ = q-p;
//        for (k = 0; k < q; k++) {
            relabun[i] = CPATT[k]; //calculates relative abundances
            fracabun[i] = CPATT[k] / sum; //calculates fractional abundances
            masslist[i] = mass + c13*k; //masses of isotopes
            //masslist[i] = mass + k; //masses of isotopes
            avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
            i = i + 1;
        }
        //masses[1] = avgMass; //avg. mass
        //masses[2] = mass + dif; //calculate mass at maximum intensity
        //masses[3] = dif - p + 1; //calculate which ion has the max intensity

        //calculates beginning mass for integration
        ii = 0;

	double[] tempArr = new double[100];
	int tempIndex=0;
	for(ii=element[6];ii>0;ii--)
	    tempArr[tempIndex++] = masslist[0]-ii*c13;
	
	for(ii=0;ii<masslist.length;ii++)
	{
	    if(masslist[ii]<=0)
		break;
	}
	     
	for(ii=0;ii<masslist.length;ii++)
	{
	    if(masslist[ii]<=0 || relabun[element[6] + ii]<=0)
		break;

	    tempArr[element[6]+ii] = masslist[ii];
	}

	masslist = tempArr;

        
        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > thres) {
                startMass = masslist[ii] - 0.5;
                break;
            }
        }
        
//        masses[4] = beginmass;
        //this.startMass = beginmass;

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > thres) {
                endMass = masslist[ii] + 0.5;
            }
        }

/*
	for(int u=0;u<masslist.length;u++)
	{
	    System.out.println(u + "mass==>>\t" + masslist[u] + "\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}
	for(int u=0;u<relabun.length;u++)
	{
	    System.out.println(u + "abun==>>\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}
	
//        masses[5] = endmass;
        this.endMass = endmass;
//        System.out.println(endmass);

        /*
System.out.println("masslist");
        for (int ii0 = 0; ii0 < masslist.length; ii0++) {
            System.out.println(masslist[ii0]);

        }


        for (int ii0 = 0; ii0 < masses.length; ii0++) {
            System.out.println("mass " + ii0 + "\t" + masses[ii0]);

        }
*/
    }

    public double getEachAbund(int i, int j)
    {
        return abund[i][j];
    }

    public double getStartMass()
    {
        return startMass;
    }

    public double getEndMass()
    {
        return endMass;
    }

    public double[] getMasslist()
    {
        return masslist;
    }
	
    public double[] getRelAbun()
    {
        return relabun;
    }

    public double[] getHighMassList()
    {
/*
	for(int u=0;u<masslist.length;u++)
	{
	    System.out.println(u + "mass==>>\t" + masslist[u]);
	    if(masslist[u]==0)
		break;
	}

	for(int u=0;u<relabun.length;u++)
	{
	    System.out.println(u + "abun==>>\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}

*/
	int start=-1, end=-1;

	for(int i=0;i<relabun.length;i++)
	{
	    if(relabun[i]>=ISODIST_THRESHOLD)
	    {
		start = i;
		break;
	    }
	}

	for(int i=start;i<relabun.length;i++)
	{
	    if(relabun[i]<=ISODIST_THRESHOLD)
	    {
		end = i;		
		break;
	    }
	}

	double[] arr = new double[end-start];

	int index=0;
	for(int i=start;i<end;i++)
	{
	    arr[index] = masslist[i];
	    index++;
	}

/*
	for(int u=0;u<arr.length;u++)
	{
	    System.out.println(u + "see==>>\t" + arr[u]);
	}
*/
	return arr;
    }
    
    public double getAvgMass() {
        return avgMass;
    }

/*
    public void setEnrichment(double enrichment)
    {
	this.enrichment = enrichment;
    }
*/

    public void calcDist(int[] element, boolean isLight)
    {
        this.element = element;

        if(isLight)
            calculateLight();
        else
            calculateHeavy();
    }
}

