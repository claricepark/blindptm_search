/**
 * @file UniProtProtein.java
 * This is the source file for edu.scripps.pms.util.spectrum.UniProtProtein
 * @author Tao Xu
 * @author Robin Park
 * @date $Date: 2010/10/29 04:16:16 $
 */



package edu.scripps.pms.util.seq;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.*;

public class UniProtProtein {
    private String output = null;
    private ArrayList<String> text = new ArrayList<String>(500);
    private ArrayList<String> acs = new ArrayList<String>(500);
    private String accession = null;
    private String sequence = null;
    private String description = "";
            // this might be problematic, some entries have more than 1 DE
    private String id = null;
     
    public UniProtProtein(String idline) {
        text.add(idline);
        id = idline;
    }
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    public String getSequence() {
        if(sequence == null) {
            populateSeq();
        }
        return sequence;
    }
    private void populateSeq() {

        boolean seqStarted = false;
        StringBuffer sb = new StringBuffer(1000);
        for(Iterator<String> it = text.iterator(); it.hasNext();) {
            String line = it.next();
            if(seqStarted) {
                line = line.trim();
                line = line.replaceAll(" ", "");
                line = line.toUpperCase();
                if(!line.startsWith("//")) {
                    sb.append(line);
                }
            } else if(line.startsWith("SQ")) {
                seqStarted = true;
            }
        }
        sequence = sb.toString();
    }
    public void setDescription(String de) {
        description = de;
    }
    public void setAccession(String ac) {
        accession = ac;
    }
    public void setId(String id) {
        this.id = id;
    }
    public ArrayList<String> getAcs() {
        return acs;
    }
    public void addLine(String line) {
        if(line.startsWith("AC")) {
            String arr [] = line.split(" +");
            for(int i = 1; i < arr.length; i++) {
                String s = arr[i];
                String temp = s.substring(0, s.length()-1);
                if(accession == null) {
                    accession = temp;
                }
                acs.add(temp);
            }
        } else if(line.startsWith("DE")) {
            // this might be problematic, some entries have more than 1 DE
            description += " " + line.split("   ")[1];
             
        }
        text.add(line);
    }
    private void formatOutput() {
        StringBuffer sb = new StringBuffer(5000);
        for(Iterator<String> it = text.iterator(); it.hasNext();) {
            sb.append(it.next());
            sb.append("\n");
        }
        
        output = sb.toString();
    }
    public String getOutput() {
        if(output == null) {
            formatOutput();
        }
        return output;
    }
    public String getSubcellularLocation() {
        StringBuffer sb = new StringBuffer(500);
        for(Iterator<String> it = text.iterator(); it.hasNext();) {
            String s = it.next();
            if(s.startsWith("CC   -!- SUBCELLULAR LOCATION")) {
                sb.append(s);
                while (it.hasNext()) {
                    String line = it.next();
                    if(line.startsWith("CC   -!-")) {
                        break;
                    } else {
                        sb.append(" " + line);
                    }
                }
            }
        }        
        return sb.toString();
    }
    public String outputFasta() {
        StringBuffer sb = new StringBuffer(300); 
        sb.append(">" + accession + "\t");
        for(Iterator<String> it = acs.iterator(); it.hasNext();) {
            sb.append(it.next()); 
            sb.append(";"); 
        }
   
        sb.append(" " + description);
        return sb.toString() + "\n" + getSequence();
    }
}
