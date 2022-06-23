package edu.scripps.pms.msfilter;

import java.io.*;
import java.util.*;

public class Cluster {

public static void normalize(int nrows, int ncolumns, double data[][])

{ int i, j;
      for (i = 0; i < nrows; i++)
      { int n = 0;
        double average = 0.0;
        for (j = 0; j < ncolumns; j++)
          { average += data[i][j];
            n++;
          }
        average /= n;
        for (j = 0; j < ncolumns; j++)
          data[i][j] -= average;
      }
}

/* ********************************************************************* */

public static
double correlation (int n, double data1[][], double data2[][], 
  int index1, int index2)
/*
Purpose
=======

The correlation routine calculates the weighted Pearson distance between two
rows or columns in a matrix. We define the Pearson distance as one minus the
Pearson correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row. 

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

index1     (input) int
Index of the first row.

index2     (input) int
Index of the second row.
============================================================================
*/
{ double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  int i;
  for (i = 0; i < n; i++) {
        double term1 = data1[index1][i];
        double term2 = data2[index2][i];
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
  }
  result -= sum1 * sum2 / n;
  denom1 -= sum1 * sum1 / n;
  denom2 -= sum2 * sum2 / n;
  if (denom1 <= 0) return 1.0; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1.0; /* include '<' to deal with roundoff errors */
  result = result / Math.sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

public static
void randomassign (int nclusters, int nelements, int clusterid[])
/*
Purpose
=======

The randomassign routine performs an initial random clustering, needed for
k-means or k-median clustering. Elements (genes or microarrays) are randomly
assigned to clusters. First, nclust elements are randomly chosen to be assigned
to the clusters 0..nclust-1 in order to guarantee that none of the clusters
are empty. The remaining elements are then randomly assigned to a cluster.

Arguments
=========

nclust  (input) int
The number of clusters.

nelements  (input) int
The number of elements to be clustered (i.e., the number of genes or microarrays
to be clustered).

clusterid  (output) int[nelements]
The cluster number to which an element was assigned.

External Functions:
ranlib: int genprm
============================================================================
*/

{ int i;
  int map[] = new int[nelements];
  /* Initialize mapping */
  for (i = 0; i < nelements; i++) map[i] = i;
  /* Create a random permutation of this mapping */
  Random generator = new Random();
  for (i=0; i<nelements; i++) {
    int randomPosition = generator.nextInt(nelements);
    int temp = map[i];
    map[i] = map[randomPosition];
    map[randomPosition] = temp;
  }

  /* Assign each of the first nclusters elements to a different cluster
   * to avoid empty clusters */
  for (i = 0; i < nclusters; i++) clusterid[map[i]] = i;

  /* Assign other elements randomly to a cluster */
  for (i = nclusters; i < nelements; i++)
    clusterid[map[i]] = generator.nextInt(nclusters);
  return;
}

/* ********************************************************************* */
public static
void getclustermean(int nclusters, int nrows, int ncolumns,
  double data[][], int clusterid[], double cdata[][])
/*
Purpose
=======

The getclustermean routine calculates the cluster centroids, given to which
cluster each element belongs. The centroid is defined as the mean over all
elements for each dimension.

Arguments
=========

nclusters  (input) int
The number of clusters.

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) double[nrows][ncolumns]
The array containing the gene expression data.

clusterid  (output) int[nrows] if transpose==0
The cluster number to which each element belongs. If transpose==0, then the
dimension of clusterid is equal to nrows (the number of genes). Otherwise, it
is equal to ncolumns (the number of microarrays).

cdata      (output) double[nclusters][ncolumns] if transpose==0
On exit of getclustermean, this array contains the cluster centroids.

========================================================================
*/
{ int i, j, k;
  int[][] cmask = new int[nclusters][ncolumns];

   for (i = 0; i < nclusters; i++)
    { for (j = 0; j < ncolumns; j++)
      { 
	cmask[i][j] = 0;
        cdata[i][j] = 0.;
      }
    }
    for (k = 0; k < nrows; k++)
    { i = clusterid[k];
      for (j = 0; j < ncolumns; j++)
      {
          cdata[i][j]+=data[k][j];
	  cmask[i][j]++;
      }
    }
    for (i = 0; i < nclusters; i++)
    { for (j = 0; j < ncolumns; j++)
        { 
	cdata[i][j] /= cmask[i][j];
        }
    }
  return;
}

/* ---------------------------------------------------------------------- */
                                                                                
public static
boolean equal_clusters(int n, int clusterids1[], int clusterids2[])
/*
This function checks if two k-means clustering solutions are equal to each
other. If equal, the function returns true; otherwise, it returns false.
                                                                                
n          (input) int
The size of the arrays clusterids1 and clusterids2, equal to the number of
items that were clustered.
                                                                                
clusterids1 (input) int[n]
An array containing n elements, indicating the number of the cluster to which
each of the items was assigned in the first clustering solution.
                                                                                
clusterids2 (input) int[n]
An array containing n elements, indicating the number of the cluster to which
each of the items was assigned in the second clustering solution.
*/
{ int i;
  for (i = 0; i < n; i++)
    if (clusterids1[i]!=clusterids2[i]) return false;
  return true;
}
                                                                                
/* ********************************************************************* */

public static
void emalg (int nclusters, int nitems, int ndata,
  double data[][], int clusterid[])

{ int cn[] = new int[nclusters];
  /* This will contain the number of elements in each cluster. This is needed
   * to check for empty clusters.
   */

  int savedids[] = new int[nitems];
  /* needed to check for periodic behavior */

  boolean changed = true;
  int iteration = 0;
  int period = 10;
  int i, j;

  for (i = 0; i < nitems; i++) cn[clusterid[i]]++;

  /* Start the loop */
  do 
  { if (iteration % period == 0)
    { /* save the current clustering solution */
      for (i = 0; i < nitems; i++) savedids[i] = clusterid[i];
      period *= 2;
    }
    iteration += 1;

    /* Find the center */
    double cdata[][] = new double[nclusters][ndata];
    Cluster.getclustermean(nclusters, nitems, ndata, data, clusterid,
                     cdata);
    for (i = 0; i < nclusters; i++) {
	for (j = 0; j < ndata; j++) {
	    data[i+nitems][j] = cdata[i][j];
	}
    }
    changed = false;

    for (i = 0; i < nitems; i++)
    /* Calculate the distances */
    { double distance;
      int jnow = clusterid[i];
      if (cn[jnow]==1) continue;
      /* No reassignment if that would lead to an empty cluster */
      /* Treat the present cluster as a special case */
      distance = Cluster.correlation(ndata,data,data,i,nitems+jnow);
      for (j = 0; j < nclusters; j++)
      { double tdistance;
        if (j==jnow) continue;
        tdistance = Cluster.correlation(ndata,data,data,i,nitems+j);
        if (tdistance < distance)
        { distance = tdistance;
          cn[clusterid[i]]--;
          clusterid[i] = j;
          cn[j]++;
          changed = true;
        }
      }
    }
  } while (changed && !Cluster.equal_clusters(nitems, savedids, clusterid));
  return;
}

/* *********************************************************************** */

public static int kcluster (int nclusters, int nrows, int ncolumns,
  double data[][], int npass, int clusterid[])
/*
Purpose
=======

The kcluster routine performs k-means or k-median clustering on a given set of
elements, using the specified distance measure. The number of clusters is given
by the user. Multiple passes are being made to find the optimal clustering
solution, each time starting from a different initial clustering.


Arguments
=========

nclusters  (input) int
The number of clusters to be found.

data       (input) double[nrows][ncolumns]
The array containing the data of the elements to be clustered (i.e., the gene
expression data).

nrows     (input) int
The number of rows in the data matrix, equal to the number of genes.

ncolumns  (input) int
The number of columns in the data matrix, equal to the number of microarrays.

npass      (input) int
The number of times clustering is performed. Clustering is performed npass
times, each time starting from a different (random) initial assignment of 
genes to clusters. The clustering solution with the lowest within-cluster sum
of distances is chosen.
If npass==0, then the clustering algorithm will be run once, where the initial
assignment of elements to clusters is taken from the clusterid array.

clusterid  (output; input) int[nrows] if transpose==0
                           int[ncolumns] if transpose==1
The cluster number to which a gene or microarray was assigned. If npass==0,
then on input clusterid contains the initial clustering assignment from which
the clustering algorithm starts. On output. it contains the clustering solution
that was found.
========================================================================
*/
{ int nelements = nrows;
  int ndata = ncolumns;

  int i, j;
  int ipass, ifound;
  int[] tclusterid;
  int[] mapping;
  double error;

  if (nelements < nclusters)
  { 
    return 0;
  }
  /* More clusters asked for than objects available */

  /* Set the result of the first pass as the initial best clustering solution */
  ifound = 1;

  double newdata[][] = new double[nelements+nclusters][ndata];
    for (i = 0; i < nelements; i++)
    { 
	for (j = 0; j < ndata; j++)
	{
	    newdata[i][j] =  data[i][j];
	}
    }
    data = newdata;
  
//  Cluster.normalize(nelements, ndata, data);

  /* Find out if the user specified an initial clustering */
  if (npass!=0)
  { /* First initialize the random number generator */
    Cluster.randomassign (nclusters, nelements, clusterid);
  }

  error = 0.;
  Cluster.emalg(nclusters, nelements, ndata, data, clusterid);

  for (i = 0; i < nelements; i++)
  { j = clusterid[i];
    error += Cluster.correlation(ndata, data, data, i, nelements+j);
  }

  tclusterid = new int[nelements];
  mapping = new int[nclusters];

  for (ipass = 1; ipass < npass; ipass++)
  { double tssin = 0.;
    boolean same = true;

    Cluster.randomassign (nclusters, nelements, tclusterid);
    Cluster.emalg(nclusters, nelements, ndata, data, tclusterid);

    for (i = 0; i < nclusters; i++) mapping[i] = -1;
    for (i = 0; i < nelements; i++)
    { j = tclusterid[i];
      if (mapping[j] == -1) mapping[j] = clusterid[i];
      else if (mapping[j] != clusterid[i]) same = false;
      tssin +=
        Cluster.correlation(ndata, data, data, i, nelements+j);
    }
    if (same) (ifound)++;
    else if (tssin < error)
    { ifound = 1;
      error = tssin;
      for (i = 0; i < nelements; i++) clusterid[i] = tclusterid[i];
    }
  }

  return ifound;
}

}
