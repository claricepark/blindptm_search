/**
 * @file Cluster3SiblingGeneRetriever.java
 * This is the source file for edu.scripps.pms.util.spectrum.Cluster3SiblingGeneRetriever
 * @author Tao Xu
 * @date $Date
 */

package edu.scripps.pms.cluster3;

import java.util.*;
import java.io.*;


// this program expect a header line in all files to be joined
public class Cluster3SiblingGeneRetriever {
    private String line; // one entry in the cluster 3.0 .gtr file
    private Cluster3SiblingGeneRetriever node1;
    private Cluster3SiblingGeneRetriever node2;
    private boolean isGene = false;
    private Cluster3SiblingGeneRetriever parent;
   
    private String nodestring;
    private String node1string;
    private String node2string;
    private String score;

    
    public Cluster3SiblingGeneRetriever(String l) {
        line = l;
        String [] arr = line.split("\t");
        nodestring = arr[0];
        isGene = nodestring.startsWith("GENE");
        node1string = arr[1];
        node2string = arr[2];
        score = arr[3];
    }
   
    public Cluster3SiblingGeneRetriever getParentCluster3SiblingGeneRetriever() {
        return parent;
    } 
    public Cluster3SiblingGeneRetriever getSiblingCluster3SiblingGeneRetriever() {
        if(parent != null) {
            if(parent.node1 != this) {
                return parent.node1;
            } else {
                return parent.node2;
            }
        } else {
            return null;
        }
    } 

    // java edu.scripps.pms.cluster3.Cluster3SiblingGeneRetriever clustered_jb_20121217.gtr clustered_jb_20121217.cdt interested_gene_symbols.txt > out.txt
    public String getSiblingGenes(String gid, HashMap<String, String> genemap) {
        if(gid.equals(node1string)) {
            if(node2string.startsWith("GENE")) {
                //return genemap.get(node2string);
                return getGeneInfo(node2string, genemap);
            } else {
                return node2.getGenes(genemap);
            }
        } else {

            if(node1string.startsWith("GENE")) {
                //return genemap.get(node1string);
                return getGeneInfo(node1string, genemap);
            } else {
                return node1.getGenes(genemap);
            }

        }

    } 
   
    // the geneid should startw with GENE, not NODE
    private String getGeneInfo(String geneid, HashMap<String,String> genemap) {
        String geneline = genemap.get(geneid);
        String genesymbol = getGeneSymbol(geneline);        
        return genesymbol + "\t" + geneid + "\t" + geneline + "\n";

    }
    public String getGenes(HashMap<String, String> genemap) { // genemap is from cluster 3.0 .cdt file
        StringBuffer sb = new StringBuffer();
        if(node1string.startsWith("GENE")) {
            //sb.append(genemap.get(node1string));
            //sb.append("\n");
            sb.append(getGeneInfo(node1string, genemap));
            if(node2string.startsWith("GENE")) {
                //sb.append(genemap.get(node2string));
                sb.append(getGeneInfo(node2string, genemap));
                //sb.append("\n");
            } else {
                sb.append(node2.getGenes(genemap));
                //sb.append("\n");
            }

        } else {
            
            sb.append(node1.getGenes(genemap));
            sb.append("\n");
            if(node2string.startsWith("GENE")) {
                //sb.append(genemap.get(node2string));
                sb.append(getGeneInfo(node2string, genemap));
                //sb.append("\n");
            } else {
                sb.append(node2.getGenes(genemap));
                //sb.append("\n");
            }

        } 

        return sb.toString();
    }
    public void assignCluster3SiblingGeneRetrievers(HashMap<String, Cluster3SiblingGeneRetriever> nodemap, HashMap<String, String> parentmap) {
        node1 = nodemap.get(node1string);
        node2 = nodemap.get(node2string);
//System.out.println("nodestring: " + nodestring);

//System.out.println("node1stirng: " + node1string + "\tnode2sring " + node2string);
//System.out.println("node1: " + node1 + "\tnode2 " + node2);

        //node1.parent = nodemap.get(parentmap.get(nodestring)); 
        if(node1 != null) {
            node1.parent = this; 
        }
        if(node2 != null) {
            node2.parent = this; 
        }
    } 

    public static String getGeneSymbol(String line) {
        String genenamestring = line.indexOf("Gene_Symbol=") > -1? "Gene_Symbol=" : "GN=";
        if(line.indexOf(genenamestring) == -1) return null; // neither Gene_Symbol= nor GN= is found.
        //String [] arr = line.split("Gene_Symbol=");
        String [] arr = line.split(genenamestring);
        if(arr.length > 1) {
            String geneele = arr[1];
            String genename = geneele.split("\\s+")[0];
            return genename;

        } else {

//System.out.println("no gene nmae: "  +  defline);
            return null;
        }
         

    }
    public static void main (String [] args) throws IOException {

        String gtrfile = args[0];
        String cdtfile = args[1];
        String targetgenefile = args[2];
 
      
        HashMap<String, String> genemap = new HashMap(); 
        HashMap<String, ArrayList<String>> genesymbol2gids = new HashMap(); 
        HashMap<String, String> gid2genesymbol = new HashMap(); 
        HashMap<String, Cluster3SiblingGeneRetriever> nodemap = new HashMap(); 
        HashMap<String, String> parentmap = new HashMap(); 
        ArrayList<Cluster3SiblingGeneRetriever> nodes = new ArrayList();

        BufferedReader br1 = new BufferedReader(new FileReader(cdtfile));

        String line = br1.readLine();   // get gene list from cdf file
        line = br1.readLine();   // get gene list from cdf file
        line = br1.readLine();   // get gene list from cdf file
        line = br1.readLine();   // get gene list from cdf file
        while (line != null) {
            if(!line.equals("")) {
                String [] arr = line.split("\t");
                genemap.put(arr[0], arr[1]);

                String genesymbol = getGeneSymbol(line);
                gid2genesymbol.put(arr[0], genesymbol);
                ArrayList<String> gids = genesymbol2gids.get(genesymbol);
                if(gids == null) {
                    gids = new ArrayList();
                    genesymbol2gids.put(genesymbol, gids);
                }
                gids.add(arr[0]);
//System.out.println(arr[0] + " added for " + genesymbol);
            }
            line = br1.readLine();
        }
       
        br1.close();
  

        BufferedReader br = new BufferedReader(new FileReader(gtrfile));
        line =  br.readLine();   // get node info from gtr file
        while (line != null) {
            if(!line.equals("")) {
//System.out.println("line in gtr file: " + line);
                Cluster3SiblingGeneRetriever n = new Cluster3SiblingGeneRetriever(line);
                nodes.add(n);
                String arr [] = line.split("\t");
                nodemap.put(arr[0], n);
                parentmap.put(arr[1], arr[0]);
                parentmap.put(arr[2], arr[0]);
                 
            }
            line = br.readLine();
        }
        br.close();

        for(Iterator<Cluster3SiblingGeneRetriever> it = nodes.iterator(); it.hasNext();) {
            Cluster3SiblingGeneRetriever n = it.next();
            n.assignCluster3SiblingGeneRetrievers(nodemap, parentmap);
        }

        BufferedReader br2 = new BufferedReader(new FileReader(targetgenefile));
        line =  br2.readLine();   // get the target protein gene symbols 
        while (line != null) {
            line = line.trim();
            if(!line.equals("")) {
                StringBuffer sb = new StringBuffer();
                System.out.println("Now processing " + line);
             
                sb.append(line);
                ArrayList<String> gids =  genesymbol2gids.get(line);
                if(gids == null) {
                    System.out.println(line + " cannot be found in the clustering result");
                } else {
                    for(Iterator<String> it = genesymbol2gids.get(line).iterator(); it.hasNext();) {
                        String gid = it.next();
                    //System.out.println("gid: " + gid);
                        sb.append("\t" + gid + "\t" + genemap.get(gid));
                        Cluster3SiblingGeneRetriever parent = nodemap.get(parentmap.get(gid));
                        sb.append("\n");
                        sb.append(parent.getSiblingGenes(gid, genemap));
                    }
                }
                sb.append("\t");
                
 
                System.out.println(sb.toString() + "\n");
         
            }
            line = br2.readLine();
        }
        br2.close();
    }
     
}
