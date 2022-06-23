/**
 * @file SearchParams.java
 * This is the source file for edu.scripps.topdown.SearchParams
 * @author Tao Xu
 * @date $Date: 2007/12/02 18:42:03 $
 */
package edu.scripps.pms.topdown;
import java.util.*;
import java.io.File;
import java.io.IOException;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
//import org.jdom.Namespace;
import org.jdom.input.SAXBuilder;
import edu.scripps.pms.util.enzyme.Protease;

public class SearchParams {
    private String databasePath;
    private String databaseName;
    private int primaryScoreType;
    private int secondaryScoreType;
    private int minMatch;
    private String precursorIsotope;
    private String fragmentIsotope;
    private int numIsotopicPeaks;
    private int preprocess;   // preprocess mode
    private double highPrecursorTolerance;
    private double lowPrecursorTolerance;
    private double fragmentTolerance;
    private double precursorTolerance;
    private double maxPrecursorMass;
    private double minPrecursorMass;
    private int minimumProteinLength;
    private int maximumProteinLength;

    private int maxNumSpectra;
    private int minNumSpectra;
    private int maxAlter;
    private boolean isDatabaseIndexed;
    private int peakRankThreshold;
    private int candidatePeptideThreshold;
//    private Modifications mods = new Modifications(); 

    private String paramFilePath;
    private String paramFileName;
 
    private Protease p;
    private int enzymeSpecificity; 
    private ArrayList<Modification> staticMods = new ArrayList<Modification>();
    private ArrayList<DiffMod> diffMods = new ArrayList<DiffMod>();
    private ArrayList<TerminalModification> ntermDiffMods = new ArrayList<TerminalModification>();
    private ArrayList<TerminalModification> ctermDiffMods = new ArrayList<TerminalModification>();
    
    private double ntermStaticMod;
    private double ctermStaticMod;
    private ArrayList<Modifications> allModifications = new ArrayList<Modifications>(1000);
    private HashSet<Modifications> mods = new HashSet<Modifications>();    
    private double maxMassShift = -10000;
    private int maxNTermTruncation;
    private int maxCTermTruncation;
    private double maxNTermMassAddition;
    private double maxCTermMassAddition;
    private double maxNTermMassLose;
    private double maxCTermMassLose;


    public SearchParams(String path, String paramFile) throws IOException, JDOMException {
        paramFilePath = path;
        paramFileName = paramFile;
        readParams();
    }

    
    public SearchParams(String paramFile) throws IOException, JDOMException {
        paramFileName = paramFile;
        paramFilePath = ".";
        readParams();
    }
    public int getMaxCTermTruncation() {
        return maxCTermTruncation;
    }
    public int getMaxNTermTruncation() {
        return maxNTermTruncation;
    }
    public double getMaxCTermMassAddition() {
        return maxCTermMassAddition;
    }
    public double getMaxNTermMassAddition() {
        return maxNTermMassAddition;
    }
    public double getMaxCTermMassLose() {
        return maxCTermMassLose;
    }
    public double getMaxNTermMassLose() {
        return maxNTermMassLose;
    }
    protected Modifications getStaticTerminalModifications() {

        Modifications temp = new Modifications();
        temp.setStaticNTermMod(ntermStaticMod);
        temp.setStaticCTermMod(ctermStaticMod);
        return temp;
    }
    
    public int getNumModifications() {
        return allModifications.size();
    }    
    public double getMaxMassShift() {
        if(maxMassShift < -9999 ) {
                 
            for(Iterator<Modifications> it = getAllModifications(); it.hasNext();) {
                double shift = it.next().getMassShift();
                maxMassShift = maxMassShift > shift? maxMassShift : shift;
            }
            maxMassShift = maxMassShift < -9999? 0 : ntermStaticMod + ctermStaticMod; // in case no modification 
        }
        return maxMassShift;

    }
    public Iterator<Modifications> getAllModifications() {
        return allModifications.iterator();
    }
    public double getStaticCTermMod() {
        return ctermStaticMod;
    }
    public double getStaticNTermMod() {
        return ntermStaticMod;
    }
    // temporary solution
    public Modifications getModifications() {
        return allModifications.get(0);
    }
    private void addAllModifications() {
        int numDiffMods = diffMods.size(); 
        Modifications m = getStaticTerminalModifications(); 
        addTerminalDiffMods(m); // with terminal diff mods only
        for(int diffModIndex = 0; diffModIndex < numDiffMods; diffModIndex++) {
            addDiffMods(diffModIndex, maxAlter, m);   
        }
        for(Iterator<Modifications> it = allModifications.iterator(); it.hasNext();) {
            System.out.println(it.next().getInfo());
        }
    }
    private void addModifications(Modifications m) {
        //System.out.print("begin in addModifications, massShift: " + m.getDiffModsShift() + "\tnumModificationsAddeSoFar: " + allModifications.size() + "\t");
        
        if(m != null && m.getDiffModsShift() != 0 && !allModifications.contains(m)) {
            allModifications.add(m.copy());
        } else {
            //System.out.println(m.getInfo() + "Not added!!!");
        } 
    }

    // add c and n terminal modifications and then add m to allModifications
    private void addTerminalDiffMods(Modifications mToAdd) {
        Modifications m = mToAdd.copy();
        int numNtermMods = ntermDiffMods.size();
        int numCtermMods = ctermDiffMods.size();
        addModifications(m); // no n no c
        if(numNtermMods == 0) {
            if(numCtermMods == 0) {
                return;
            } else {
                for(int i = 0; i < numCtermMods; i++) {
                    Modifications newm = m.copy();
                    newm.setCTerm(ctermDiffMods.get(i));
                    addModifications(newm);
                }
            }
        } else {
            if(numCtermMods == 0) {
                for(int i = 0; i < numNtermMods; i++) {
                    Modifications newm = m.copy();
                    newm.setNTerm(ntermDiffMods.get(i)); 
                    addModifications(newm);
                }
            } else {
                for(int i = 0; i < numNtermMods; i++) {
                    Modifications newm = m.copy();
                    for(int j = 0; j < numCtermMods; j++) {
                        Modifications nc = m.copy(); 
                        nc.setCTerm(ctermDiffMods.get(j));
                        addModifications(nc);   // no n + c
                    }                     
                    addModifications(newm);
                    newm.setNTerm(ntermDiffMods.get(i)); // n no c
                    for(int j = 0; j < numCtermMods; j++) {
                        Modifications nc = newm.copy();
                        addModifications(nc);       // n + c
                        nc.setCTerm(ctermDiffMods.get(j));
                    }
                }
            }
        } 
        //System.out.println("end in addTerminalDiffMods");
    }
    private void addDiffMods(int index, int numToAdd, Modifications temp) {
        if(numToAdd < 1 || index >= diffMods.size() || temp.getNumDiffMods() >= maxAlter) {
            return;
        }
        //addModifications(m); 
        Modifications newm = temp.copy();
        //System.out.println("in addDiffMods, indiex: " + index + "\tnumToAdd: " + numToAdd + "\tnumModsInTemp: " + temp.getInfo());
        //addTerminalDiffMods(m); // this will add  m and its siblings ith terminal 


        addDiffMods(index+1, numToAdd, newm); // without add this one
        for(int i = 1; i <= numToAdd; i++) {
            newm = newm.copy();
            newm.addDiffMod(diffMods.get(index)); // add one
            addTerminalDiffMods(newm);
            addDiffMods(index+1, numToAdd-i, newm); // after add one

        }
    }
    private void readParams() throws IOException, JDOMException {
        File paraFile = new File(paramFilePath + "/" + paramFileName); 
        Document doc = new SAXBuilder().build(paraFile);
        
        Element root = doc.getRootElement();
        readDatabase(root.getChild("database"));
        readSearchMode(root.getChild("search_mode"));
        readIsotopes(root.getChild("isotopes"));
        readTolerance(root.getChild("tolerance"));
        readPrecursorMassLimits(root.getChild("precursor_mass_limits"));
        readProteinLengthLimits(root.getChild("protein_length_limits"));
        readNumSpectraLimits(root.getChild("num_spectra_limits"));
        maxAlter =  Integer.parseInt(root.getChildTextTrim("max_alter")); 
        readModifications(root.getChild("modifications"));
        readEnzymeInfo(root.getChild("enzyme_info"));

    }
    private void readEnzymeInfo(Element e) {
        enzymeSpecificity = Integer.parseInt(e.getChildTextTrim("specificity")); 
        if(enzymeSpecificity > 0) {
            String enzymeName = e.getChildTextTrim("name");
            p = new Protease(enzymeName);
            p.setType(Boolean.parseBoolean(e.getChildTextTrim("type")));
            Element re = e.getChild("residues"); 
            List residues = re.getChildren();
            for(Iterator i = residues.iterator(); i.hasNext();) {
                Element r = (Element)i.next();
                String s = r.getTextTrim();
                p.addCleavageSite(s.charAt(0));
            }
        }
    }
    private void readDatabase(Element e) {
        databaseName = e.getChildTextTrim("database_name");
        //databasePath = e.getChildTextTrim("database_path");
        //if(databasePath != null) {
           // databaseName = (databasePath + "/" + databaseName).trim();
        //}
//System.out.println("\n\n\n\ndatabase_name: " + databaseName);
        String s = e.getChildTextTrim("is_indexed");
        if(s != null) {
            char c = s.charAt(0);
            if(c == 'T' || c == 't' || c == 'Y' || c == 'y') {
                isDatabaseIndexed = true;
            } else if(c == 'N' || c == 'n' || c == 'F' || c == 'f') {
                isDatabaseIndexed = false;
            } else {
                throw new InvalidArgumentException(
                    "is_indexed unknown. Should be true or false, or yes or no.");
            } 
        } else {
            isDatabaseIndexed = false; // by default, database is not indexed
        }
        //isDatabaseIndexed = Boolean.parseBoolean(e.getChildTextTrim("is_indexed"));
    }
    private void readTerminalMods(Element me, String terminal) {
//System.out.println(terminal);
        boolean isCterm = terminal.startsWith("c");
        Element e = me.getChild(terminal);
        if(e != null) {
            Element staticMod = e.getChild("static_mod");
            if(staticMod != null) {
//System.out.println("in readTerminalMods, terminal: " + terminal + "\t" + staticMod.getChildTextTrim("mass_shift"));
                double massShift = Double.parseDouble(staticMod.getChildTextTrim("mass_shift")); 
                if(isCterm) {
                    ctermStaticMod = massShift; 
                } else {
                    ntermStaticMod = massShift; 
                }
            //    char symbol = e.getChildTextTrim("symbol").charAt(0);
            }
//System.out.println("in readTerminalMods, ctermStaticMod: " + ctermStaticMod + "\tntermStaticMod: " + ntermStaticMod);
            Element diffMods = e.getChild("diff_mods");
            if(diffMods != null) {
                List <Element> diffModList = diffMods.getChildren();
                for(Element de : diffModList) {
                    //char residue = e.getChildTextTrim("residue").charAt(0);
                    char symbol = de.getChildTextTrim("symbol").charAt(0);
                    String massShift = de.getChildTextTrim("mass_shift");
                    //System.out.println("static modification: " + staticModification);
                    if (massShift != null) {
                        //m.setStaticModification(Double.parseDouble(staticModification));
                        double mass = Double.parseDouble(massShift);
                        if(mass != 0) {
                            TerminalModification d = new TerminalModification(symbol, mass);
                            if(isCterm) { 
                                ctermDiffMods.add(d);
                            } else {
                                ntermDiffMods.add(d);
                            }
                        }
                        System.out.println(terminal + " diff mod: " + massShift);
                    }
                } 

            }
        }
    }
    private void readModifications(Element me) {
       // System.out.println("In readModifications"); 
     
        maxNTermTruncation = Integer.parseInt(me.getChildTextTrim("max_nterm_truncation"));
        maxCTermTruncation = Integer.parseInt(me.getChildTextTrim("max_cterm_truncation"));
        maxNTermMassAddition = Double.parseDouble(me.getChildTextTrim("max_nterm_mass_addition"));
        maxCTermMassAddition  = Double.parseDouble(me.getChildTextTrim("max_cterm_mass_addition"));
        maxNTermMassLose  = Double.parseDouble(me.getChildTextTrim("max_nterm_mass_lose"));
        maxCTermMassLose  = Double.parseDouble(me.getChildTextTrim("max_cterm_mass_lose"));
        System.out.println("Maximum num residues for nterm truncation: " + maxNTermTruncation);
        System.out.println("Maximum num residues for cterm truncation: " + maxCTermTruncation);
        System.out.println("Maximum num residues for nterm mass addition: " + maxNTermMassAddition);
        System.out.println("Maximum num residues for cterm mass addition: " + maxCTermMassAddition);
        System.out.println("Maximum num residues for nterm mass lose: " + maxNTermMassLose);
        System.out.println("Maximum num residues for cterm mass lose: " + maxCTermMassLose);
        
        readTerminalMods(me, "c_term");
        readTerminalMods(me, "n_term");
        // get static mods
        List <Element> staticModList = me.getChild("static_mods").getChildren();
        for(Element e : staticModList) {
            char residue = e.getChildTextTrim("residue").charAt(0);
            String massShift = e.getChildTextTrim("mass_shift");
            //System.out.println("static modification: " + staticModification);
            if (massShift != null) {
                double mass = Double.parseDouble(massShift);
                if(mass != 0) {
                    staticMods.add(new Modification(residue, mass));
                }
                System.out.println("static modification: " + massShift);
            }
        } 
        // get diff_mods
        List <Element> diffModList = me.getChild("diff_mods").getChildren();
        for(Element e : diffModList) {
            //char residue = e.getChildTextTrim("residue").charAt(0);
            char symbol = e.getChildTextTrim("symbol").charAt(0);
            String massShift = e.getChildTextTrim("mass_shift");
            
            //System.out.println("static modification: " + staticModification);
            if (massShift != null) {
                //m.setStaticModification(Double.parseDouble(staticModification));
                double mass = Double.parseDouble(massShift);
                if(mass != 0) {
                    String neutralLoss = e.getChildTextTrim("neutral_loss");         
                    double nl = 0; 
                    if(neutralLoss != null) {
                        nl = Double.parseDouble(neutralLoss);
                    }
                    DiffMod d = new DiffMod(mass, nl, symbol);
                    List<Element> residues = e.getChild("residues").getChildren();
                    for(Element r : residues) {
                        byte c = (byte)r.getTextTrim().charAt(0);
                        System.out.print("residue " + (char)c + "\t mass shift " + mass + "\n");
                        d.setModifiable(c, true);
                    }
                    diffMods.add(d);
                }
                System.out.println("diff mod: " + massShift);
            }
        } 
        addAllModifications();  
        System.out.println("Number of Modifications: " + getNumModifications()); 
        System.out.println("Number of Diff Mods: " + diffMods.size()); 
    }
    private void readNumSpectraLimits(Element ne) {
        maxNumSpectra =  Integer.parseInt(ne.getChildTextTrim("maximum"));
        minNumSpectra =  Integer.parseInt(ne.getChildTextTrim("minimum"));
    }
    private void readPrecursorMassLimits(Element pe) {
//        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
//        minPrecursorMass = Double.parseDouble(pe.getChildTextTrim("minimum"));
        minPrecursorMass = Double.parseDouble(pe.getChildTextTrim("minimum"));
    }

    private void readProteinLengthLimits(Element pe) {
        minimumProteinLength = Integer.parseInt(pe.getChildTextTrim("minimum"));
        maximumProteinLength = Integer.parseInt(pe.getChildTextTrim("maximum"));
    }
 
    private void readTolerance(Element te) {
//        precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        //precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        highPrecursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor_high"));
        lowPrecursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor_low"));
//        fragmentTolerance = Double.parseDouble(te.getChildTextTrim("fragment"));
        precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        fragmentTolerance = Double.parseDouble(te.getChildTextTrim("fragment"));

    }
    private void readIsotopes(Element ie) 
                              throws InvalidArgumentException {
        String mono = MassSpecConstants.MONOISOTOPE;
        String avg = MassSpecConstants.AVGISOTOPE;
        String isotope =  ie.getChildTextTrim("precursor");
        numIsotopicPeaks = Integer.parseInt(ie.getChildTextTrim("num_peaks")); 
        // set the precursor isotope type
        if (mono.equals(isotope) || avg.equals(isotope)) {
            precursorIsotope = isotope;
        } else {
            throw new InvalidArgumentException(
                "Unknown precursor isotope type. Should be mono or avg.");
            
        }

        // set fgragment isotope type
        isotope = ie.getChildTextTrim("fragment");     
        if (mono.equals(isotope) || avg.equals(isotope)) {
            fragmentIsotope =  isotope;
        } else {
            throw new InvalidArgumentException(
                "Unknown fragment isotope type. Should be mono or avg.");
        }        
    }
    private void readSearchMode(Element me) {
        primaryScoreType = Integer.parseInt(me.getChildTextTrim("primary_score_type"));
        peakRankThreshold = Integer.parseInt(me.getChildTextTrim("peak_rank_threshold")); 
        candidatePeptideThreshold = Integer.parseInt(me.getChildTextTrim("candidate_peptide_threshold")); 
        secondaryScoreType =  Integer.parseInt(me.getChildTextTrim("secondary_score_type"));
        minMatch =  Integer.parseInt(me.getChildTextTrim("min_match"));
        preprocess  =  Integer.parseInt(me.getChildTextTrim("preprocess"));
        String fragMethod = me.getChildTextTrim("fragmentation_mothod");
        
    }
    public void setDbName(String dbName) {
        databaseName = dbName;
    }
    public void setPrimaryScoreType(int primaryScoreType) {
        this.primaryScoreType = primaryScoreType;
    }

    public void setsecondaryScoreType(int secondaryScoreType) { 
        this.secondaryScoreType = secondaryScoreType;
    }
    public void setMinMatch(int minMatch) { 
        this.minMatch = minMatch;
    }

    public void setPrecursorIsotope(String precursorIsotope) 
                              throws InvalidArgumentException {

        if (MassSpecConstants.MONOISOTOPE.equals(precursorIsotope) ||
                MassSpecConstants.AVGISOTOPE.equals(precursorIsotope)) {
            this.precursorIsotope = precursorIsotope;
        } else {
            throw new InvalidArgumentException(
                "Unknow precursor isotope type. Should be mono or avg.");
        }
    }
    public void setFragmentIsotope(String fragmentIsotope) 
                              throws InvalidArgumentException {

        if (MassSpecConstants.MONOISOTOPE.equals(fragmentIsotope) ||
                MassSpecConstants.AVGISOTOPE.equals(fragmentIsotope)) {
            this.fragmentIsotope = fragmentIsotope;
        } else {
            throw new InvalidArgumentException(
                "Unknow fragment isotope type. Should be mono or avg.");
        }
    }
   
    public void setPreprocess(int mode) {
        preprocess = mode;
    }
    public void setFragmentTolerance(double tolerance) {
        fragmentTolerance = tolerance;
    }
    public void setMaxPrecursorMass(double mass) {
        maxPrecursorMass = mass;
    }
    public void setMinPrecursorMass(double mass) {
        minPrecursorMass = mass;
    }
    public void setMinimumProteinLength(int len) {
        minimumProteinLength = len;
    }
    public void setMaxNumSpectra(int numSpectra) {
        maxNumSpectra = numSpectra;
    }
    public void setMinNumSpectra(int numSpectra) {
        minNumSpectra = numSpectra;
    }
    public void setMaxAlter(int maxAlter) {
        this.maxAlter = maxAlter;
    }
    public int getPeakRankThreshold() {
        return peakRankThreshold;
    }
    public int getCandidatePeptideThreshold() {
        return candidatePeptideThreshold;
    }
    public Iterator<Modification> getStaticMods() {
        return staticMods.iterator();
    }
    public Iterator<DiffMod> getDiffMods() {
        return diffMods.iterator(); 
    }
    public String getDbName() {
        return databaseName;
    }    
    public boolean isDatabaseIndexed() {
        return isDatabaseIndexed;
    }

    public int getPrimaryScoreType() {
        return primaryScoreType;
    }

    public int getSecondaryScoreType() {
        return secondaryScoreType;
    }
    
    public TerminalModification getNTerminalModification() {
        return ctermDiffMods.get(0);
    }
    public TerminalModification getCTerminalModification() {
        return ctermDiffMods.get(0); 
    }
    public int getMinMatch() {
        return minMatch;
    }
    public int getNumIsotopicPeaks() {
        return numIsotopicPeaks;
    }
    public String getPrecursorIsotope() {
        return precursorIsotope;
    }

    public String getFragmentIsotope() {
        return fragmentIsotope;
    }
    public int getPreprocess() {
        return preprocess;
    }

    public double getHighPrecursorTolerance() {
        return highPrecursorTolerance;
    }
    public double getLowPrecursorTolerance() {
        return lowPrecursorTolerance;
    }
    public double getPrecursorTolerance() {
        return precursorTolerance;
    }
    public double getFragmentTolerance() {
        return fragmentTolerance;
    }

    public int getMinimumProteinLength() {
        return minimumProteinLength;
    }
    public int getMaximumProteinLength() {
        return maximumProteinLength;
    }
    public double getMaxPrecursorMass() {
        return maxPrecursorMass;
    }
    public double getMinPrecursorMass() {
        return minPrecursorMass;
    }

    public int getMaxNumSpectra() {
        return maxNumSpectra;
    }
    public int getMinNumSpectra() {
        return minNumSpectra;
    }
    public int getMaxAlter() {
        return maxAlter;
    }
    public Protease getProtease() {
        return p;
    }
    public int getEnzymeSpecificity() {
        return enzymeSpecificity;
    }
    public String getEnzymeName() {
        return p.getName();
    }
    public double getStaticTerminalMods() {
        //return mods.getStaticTerminalMods();
        return ctermStaticMod + ntermStaticMod;
    }
    public static void main(String args[]) throws Exception {
        SearchParams sp = new SearchParams(args[0], args[1]);
        System.out.println("Database is: " + sp.getDbName());
        System.out.println("maxAlter is: " + sp.getMaxAlter());
        System.out.println("highPrecursorTolerance is: " + sp.getHighPrecursorTolerance());
        System.out.println("Fragment isotope is: " + sp.getFragmentIsotope());
        Iterator <Modification> modifications = sp.getStaticMods();
        while(modifications.hasNext()) {
            Modification m = modifications.next();
            System.out.println("Residue " + (char)m.getResidue() + "\t" + m.getMassShift() + "\n");
        }

        String a = "java";
        String b = new StringBuffer(a).toString();
        a = b;
        if(a.equals(b)) {
            System.out.println("a equal b");
 
        } 
        if(a == b) {

            System.out.println("a == b");
        }else {

            System.out.println("a != b");
        }
        
    }
}



