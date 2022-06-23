

/**
 * @author Tao Xu 
 * @version $Id: Pmsdb.java,v 1.4 2009/11/19 23:36:53 taoxu Exp $
 * @date $Date: 2009/11/19 23:36:53 $
 */


package edu.scripps.pms.pmsdb; 

import ubic.db.*;
import edu.scripps.pms.models.db.*;
import java.sql.*;
import java.util.Date;
import java.util.HashSet;
import java.io.*;


import org.apache.log4j.*;

/**
 * This class is a subclass of MysqlDb. It provides an interface to
 * access Pmsdb database.
 *
 */
public class Pmsdb extends MysqlDb {
    private static final String USAGE = "java Pmsdb outputFileName"; 
    public static final String HEADER = "<?xml version=\"1.0\"?>\n" +
                                  "<!DOCTYPE Experiment_List SYSTEM \"hspms.dtd\">\n" +
                                  "<Experiment_List>\n "; 
    public static final String TAIL = "</Experiment_List>\n";
    public static final String LOG = "SalivaDataExchange.log";
    private PrintWriter log; 
    private PrintWriter xmlFile;
    private String xmlFileName;
    private String paramFileName;
    private HashSet<String> proteins = new HashSet<String>(5000);
//    private HashSet<String> tempProteins = new HashSet<String>(5000);
    /**
     * @fn Pmsdb(String dbName)
     * Constucting a Pmsdb object
     * @param dbName - the name of the database to be accessed
     * @exception SQLException
     */
    public Pmsdb(String dbName) {
        super(dbName);
    }
    public Pmsdb(String dbName, String outfile) throws SQLException, IOException {
        super(dbName);
        log = new PrintWriter(new FileOutputStream(LOG, true));
        xmlFileName = outfile;
        xmlFile = new PrintWriter(outfile);

    }
    protected void printLog(String msg) {
        log.println(msg);
    } 
    protected void closeLog() {
        log.close();
    }
    protected int [] projectId2ExperimentIds(int id) throws SQLException {
   
        String query = "SELECT " +
                            "id " +
                        "FROM " +
                            "msp_experiment " +
                        "WHERE " +
                            "project_id = ?";

        PreparedStatement ps = prepareStatement(query);
        ps.setInt(1, id);

        ResultSet rs = ps.executeQuery();
        int numRows = getNumRows(rs);
        int [] experimentIds = new int[numRows];
                
        while (rs.next()) {
            experimentIds[--numRows] = rs.getInt("id");
        }
        return experimentIds; 
    }    
    // get all experiment ids that is greater than minExperimentId and belong to the project identified by pid 
    protected int [] projectId2ExperimentIds(int pid, int minExperimentId) throws SQLException {
   
        String query = "SELECT " +
                            "m.id " +
                        "FROM " +
                            "msp_experiment m, mass_spec_instrument msi " +
                        "WHERE " +
                            "project_id = ? AND m.id >= ? AND m.mass_spec_instrument_id = msi.id AND " +
                            "msi.instrument_type_id IN (SELECT it.id FROM instrument_type it)";

        PreparedStatement ps = prepareStatement(query);
        ps.setInt(1, pid);
        ps.setInt(2, minExperimentId);

        ResultSet rs = ps.executeQuery();
        int numRows = getNumRows(rs);
        int [] experimentIds = new int[numRows];
                
        while (rs.next()) {
            experimentIds[--numRows] = rs.getInt("m.id");
        }
        return experimentIds; 
    }    
    public String accession2Defline(String ac) throws SQLException {
        String query = "SELECT description FROM protein WHERE accession = '" + ac + "'";
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            return rs.getString("description");
        }
        return "";
    }
    // test if an experiment has DTASelect result or not
    protected boolean experimentHasResults(int experimentId) throws SQLException {
        String query = "SELECT id FROM db_search where msp_experiment_id = " + experimentId; 
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            if(dbSearchHasResults(rs.getInt("id"))) {
                return true;
            }
        }
        return false;
    } 

    protected boolean dbSearchHasResults(int dbSearchId) throws SQLException {
        String query = "SELECT id FROM result_selection WHERE db_search_id = " + dbSearchId;
        if (executeQuery(query).next()) {
            return true;
        }
        return false;
    }
    protected int [] experimentId2ProteinHitIds(int id) throws SQLException {
        String query = "SELECT p.id FROM protein_hit p, db_search d, result_selection r " +
                       "WHERE p.result_selection_id = r.id AND r.db_search_id = d.id AND " + 
                       "d.msp_experiment_id = " + id;
        ResultSet rs = executeQuery(query);
        int numRows = getNumRows(rs);
        int [] result = new int[numRows];

        while(rs.next()) {
            result[--numRows] = rs.getInt("p.id");
        }
        return result;
    }

    protected String experimentId2SampleSourceId(int id) throws SQLException {
        String query = "SELECT sample_type FROM sample s, msp_experiment m WHERE " + 
                       "s.id = m.sample_id AND m.id = " + id;
        ResultSet rs = executeQuery(query);
        if(rs.next()) {
            String sampleType = rs.getString("sample_type");
            if(sampleType != null) {
                sampleType = sampleType.trim();
                if(sampleType.equals("parotid")) {
                    return 2 + "";
                } else if(sampleType.equals("SM/SL")) {
                    return 5 + "";
                } else if(sampleType.equals("SM")) {
                    return 3 + "";
                } else if(sampleType.equals("SL")) {
                    return 4 + "";
                } else {
                    return 1 + "";
                }
            }
        }
        return null;

    }
    protected String experimentId2InstrumentType(int id) throws SQLException {
        String q = "SELECT instrument_type FROM instrument_type i, mass_spec_instrument m, " +
                   "msp_experiment e WHERE e.id = " + id + " AND " + 
                   "e.mass_spec_instrument_id = m.id AND m.instrument_type_id = i.id";
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getString("instrument_type");
        }
        return null;
    }

    private String sequestEnzymeNo2Enzyme(int num) {
        String e = null;
        switch(num) { 
            case 0 : e = "No_Enzyme"; break;
            case 1 : e = "Trypsin"; break;
            case 2 : e = "Chymotrypsin"; break;
            case 3 : e = "Clostripain"; break;
            case 4 : e = "Cyanogen_Bromide"; break;
            case 5 : e = "IodosoBenzoate"; break;
            case 6 : e = "Proline_Endopept"; break;
            case 7 : e = "Staph_Protease"; break; 
            case 8 : e = "Trypsin_K"; break;
            case 9 : e = "Trypsin_R"; break;
            case 10 : e = "AspN"; break;
            case 11 : e = "Cymotryp/Modified"; break;
            case 12 : e = "Elastase"; break;
            case 13 : e = "Elastase/Tryp/Chymo"; break;

        }
        return e;
    }
    protected String proteinDbId2ProteinDbName(int id) throws SQLException {
        String q = "SELECT file_name FROM protein_db WHERE id = " + id;
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getString("file_name");
        }
        return null;
    }
    protected int experimentId2ResultSelectionId(int expId) throws SQLException {
        int dbSearchId = experimentId2DbSearchId(expId);
        String q = "SELECT id FROM result_selection WHERE db_search_id = " + dbSearchId;
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getInt("id");
        }
        return -1;
    }
    protected int experimentId2DbSearchId(int id) throws SQLException {
        String q = "SELECT id FROM db_search WHERE msp_experiment_id = " + id;
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getInt("id");
        }
        return -1;
    }
    protected int experimentId2DbSearchParametersId(int id) throws SQLException {
        String q = "SELECT db_search_parameters_id as id FROM db_search  WHERE msp_experiment_id = " + id;
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getInt("id");
        }
        return -1;
    }
    protected String experimentId2DbName(int id) throws SQLException {
        int paramId = experimentId2DbSearchParametersId(id);
        String q = "SELECT file_name from protein_db p, db_search_parameters dsp " +
                   "WHERE dsp.id = " + paramId + " AND p.id = dsp.protein_db_id";
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getString("file_name");
        }
        return null;
    }
    public String outputContent(String indent, String tag, String content) {
        return indent + "<" + tag + ">" + content + "</" + tag + ">\n";
    }
    private String getMassType(int type) {
        switch(type) {
            case 0 : return "avg";
            case 1 : return "mono";
            default : return "unknown"; 
        }
    }
    private void outputDataFiles(PrintWriter pw, int experimentId, String indent) throws SQLException {
        indent += "\t";
        String q = "SELECT * FROM data_source WHERE msp_experiment_id = " + experimentId;
        ResultSet rs = executeQuery(q);
        StringBuffer sb = new StringBuffer(500);
        sb.append(indent + "<Data_Files>\n");
        String subIndent = indent + "\t";
        String moreIndent = subIndent + "\t";
System.out.println("experimentId " + experimentId);
        //String extractor = null;
        while(rs.next()) {
            sb.append(subIndent + "<Data_File>\n");            
            String ms2File = rs.getString("extracted_file_name");
            String rawFile = ms2File.substring(0, ms2File.indexOf(".ms2")) + ".RAW";
            sb.append(outputContent(moreIndent, "Data_File_ID", rs.getString("id")));
            sb.append(outputContent(moreIndent, "Raw_File_Name", rawFile)); 
            sb.append(outputContent(moreIndent, "Peaklists_File_Name", ms2File));
            String extractor = rs.getString("extractor");
            sb.append(outputContent(moreIndent, "Peaklist_Creation_Program", extractor));
            sb.append(subIndent + "</Data_File>\n");            
        }
        sb.append(indent + "</Data_Files>\n");
        pw.print(sb.toString());
    }
   
    private void msmsFractionId2Spectrum(StringBuffer sb, String id, String indent) throws SQLException {
        String msmsFractionBtag = indent + "<MS_MS_Fraction>\n";    
        String msmsFractionEtag = indent + "</MS_MS_Fraction>\n";    
        indent += "\t";
        String q = "SELECT data_source_id as did, m.m_to_z as mz, charge_state as z, loscan as s FROM " + 
                   " data_source ds, peak_list pl, msmsfraction m WHERE " + 
                   "m.id = " + id + " AND m.peak_list_id = pl.id"; 
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            sb.append(msmsFractionBtag);
            sb.append(outputContent(indent, "M_to_Z", rs.getString("mz")));
            sb.append(outputContent(indent, "Charge", rs.getString("z")));
            sb.append(outputContent(indent, "Data_File_ID", rs.getString("did")));
            sb.append(outputContent(indent, "Scan_no", rs.getString("s")));
            sb.append(msmsFractionEtag);
        }
    }
    private String peptideHitId2Ptms(String peptideHitId, String indent) throws SQLException {
        StringBuffer sb = new StringBuffer(300);
        sb.append(indent + "<PTM>\n");
        boolean hasResult = false;
        String q = "SELECT position, modified_residue, mass_shift FROM ptm p, " + 
                   "modification_type m WHERE p.peptide_hit_id = " + peptideHitId + 
                   " AND modification_type_id = m.id";
        ResultSet rs = executeQuery(q);
        String moreIndent = indent + "\t";
        while(rs.next()) {
//System.out.println("PTM found, peptideHitId: " + peptideHitId);
            hasResult = true;
            sb.append(outputContent(moreIndent, "Modified_Residue", rs.getString("modified_residue")));    
            sb.append(outputContent(moreIndent, "Position", rs.getString("position")));    
            sb.append(outputContent(moreIndent, "Mass_Shift", rs.getString("mass_shift")));    
        } 
        
        if(!hasResult) {
            return "";
        } else {
            sb.append(indent + "</PTM>\n");
            return sb.toString();
        }
    }
    private void proteinHitId2PeptideHits(StringBuffer sb, int proteinHitId, String indent) throws SQLException {

        String peptideHitsBTag = indent + "<Peptide_Hits>\n"; 
        String peptideHitsETag = indent + "</Peptide_Hits>\n"; 
       
        sb.append(peptideHitsBTag);
        indent += "\t";
        String peptideHitBTag = indent + "<Peptide_Hit>\n"; 
        String peptideHitETag = indent + "</Peptide_Hit>\n"; 
        indent += "\t";
        String q = "SELECT p.*, score, rank from peptide_hit p, dtaselect_best_peptide d, score s " + 
                   "WHERE d.protein_hit_id = " + proteinHitId +
                   " AND d.peptide_hit_id = p.id AND s.peptide_hit_id = p.id AND s.score_type_id = 1";
        ResultSet rs = executeQuery(q);
       
        while(rs.next()) {
            sb.append(peptideHitBTag);
            sb.append(outputContent(indent, "Peptide_Sequence", rs.getString("sequence"))); 
            sb.append(outputContent(indent, "Peptidehit_Score", rs.getString("primary_score")));
            sb.append(outputContent(indent, "XC", rs.getString("primary_score")));
            sb.append(outputContent(indent, "DeltaCn", rs.getString("delta_cn")));
            String peptideHitId = rs.getString("id");


                sb.append(outputContent(indent, "Sp", rs.getString("score")));
                sb.append(outputContent(indent, "Rsp", rs.getString("rank")));
/*
            String q2 = "SELECT score, rank FROM score WHERE peptide_hit_id = " + peptideHitId +  
                        " AND score_type_id = 1";
            ResultSet r = executeQuery(q2);
            if(r.next()) {
                sb.append(outputContent(indent, "Sp", r.getString("score")));
                sb.append(outputContent(indent, "Rsp", r.getString("rank")));
            }
*/
            msmsFractionId2Spectrum(sb, rs.getString("msmsfraction_id"), indent); 
            sb.append(peptideHitId2Ptms(peptideHitId, indent)); 
            sb.append(peptideHitETag);
        }
        
        sb.append(peptideHitsETag);
    }
    private String proteinHitId2ProteinInfo(int phId, String indent, String dbName) throws SQLException {
        indent += "\t";        
        String beginTag = indent + "<Protein_Hit>\n";
        String endTag = indent + "</Protein_Hit>\n";
        String q = "SELECT sequence_coverage as sc, sequence, accession, description FROM " + 
                   "protein_hit_identified_protein ph, protein p WHERE " +
                   "ph.protein_hit_id = " + phId + " AND ph.protein_id = p.id"; 
        ResultSet rs = executeQuery(q);
        String subIndent = indent + "\t";
        String moreIndent = subIndent + "\t";
        
        String proteinBTag = subIndent + "<Protein>\n"; 
        String proteinETag = subIndent + "</Protein>\n"; 
        StringBuffer sb = new StringBuffer(10000);
        while(rs.next()) { // one proteinHitId may be associated with multiple proteins (redundent)
            String acc = rs.getString("accession");
            if(acc.startsWith("Reverse")) {
                return "";
            }
            //tempProteins.add(acc);
            proteins.add(acc);
            sb.append(beginTag);

            sb.append(outputContent(subIndent, "Proteinhit_Sequence_Coverage", rs.getString("sc")));

            sb.append(proteinBTag);
            sb.append(outputContent(moreIndent, "Accession_Number", acc));
            sb.append(outputContent(moreIndent, "Protein_Sequence", rs.getString("sequence")));
            sb.append(outputContent(moreIndent, "Description", rs.getString("description").substring(1)));
            sb.append(outputContent(moreIndent, "Protein_Database_ID", dbName));
            sb.append(proteinETag);

            proteinHitId2PeptideHits(sb, phId, subIndent);
            sb.append(endTag);
            
        }
        return sb.toString();
    } 
    private int outputProteinHits(PrintWriter pw, int experimentId, String indent) throws SQLException {
        //indent += "\t";
        //tempProteins = new HashSet<String>(5000);
        String dbName = experimentId2DbName(experimentId);
        int resultSelectionId = experimentId2ResultSelectionId(experimentId);
        String q = "SELECT id FROM protein_hit WHERE result_selection_id = " + resultSelectionId;
        ResultSet rs = executeQuery(q);
        //StringBuffer sb = new StringBuffer(500);
        
        pw.print(indent + "<Protein_Hits>\n");
        int numProteins = 0; 
        while(rs.next()) {
            int proteinHitId = rs.getInt("id");
            String proteinHitInfo = proteinHitId2ProteinInfo(proteinHitId, indent, dbName);
            if(!proteinHitInfo.equals("")) {
                pw.print(proteinHitInfo);            
                numProteins++;
            }
            
        } 

        pw.print(indent + "</Protein_Hits>\n");
        //pw.print(sb.toString());
        return numProteins;
    }
    protected String dbsearchParamId2DbSearchProgram(int id) throws SQLException {

        String q = "SELECT p.* FROM program p, db_search_parameters d WHERE d.id = " + id +
                   " AND p.id = d.program_id";
        ResultSet rs = executeQuery(q);
        if(rs.next()) {
            return rs.getString("name");
        }
        return null;
    }
    private int outputDbsearch(PrintWriter pw, int experimentId, String indent) throws SQLException, IOException {
        indent += "\t";
        // output <DB_Search_Parameter_File> here
        int dbSearchParamId = experimentId2DbSearchParametersId(experimentId);
        String q = "SELECT * FROM db_search_parameters WHERE id = " + dbSearchParamId;
        ResultSet rs = executeQuery(q);

        if(rs.next()) {
            String proteindb = proteinDbId2ProteinDbName(rs.getInt("protein_db_id"));
            paramFileName = rs.getString("parameters_file");
            pw.print(outputContent(indent, "DB_Search_Parameter_File", paramFileName)); 
            // escape the & character in the parameter file
            String parameters = rs.getString("parameters"); 
            
            //pw.print(outputContent(indent, "DB_Search_Parameters", parameters)); 
            pw.print(indent + "<" + "DB_Search_Parameters" + ">\n");
            pw.print(parameters);
            pw.print(indent + "<" + "/DB_Search_Parameters" + ">\n");
            
            pw.print(outputContent(indent, "Peptide_Mass_Tolerance", rs.getString("peptide_mass_tolerance")));
            pw.print(outputContent(indent, "Fragment_Ion_Tolerance", rs.getString("fragment_ion_tolerance")));
            pw.print(outputContent(indent, "Peptide_Mass_Value_Type", getMassType(rs.getInt("mass_type_parent"))));
            pw.print(outputContent(indent, "Fragment_Mass_Value_Type", getMassType(rs.getInt("mass_type_fragment"))));
            pw.print(outputContent(indent, "DB_Search_Program_ID", dbsearchParamId2DbSearchProgram(dbSearchParamId)));
            pw.print(outputContent(indent, "Protein_Database_ID", proteindb));
            int enzymeNo = rs.getInt("enzyme_number");
            pw.print(outputContent(indent, "Enzyme", sequestEnzymeNo2Enzyme(enzymeNo)));
            pw.print(indent + "<!--Max_Missed_Cleavages will not be outputted if No_Enzyme-->\n");
            if (enzymeNo != 0) {            
                pw.print(outputContent(indent, "Max_Missed_Cleavages", rs.getString("max_missed_cleavages")));
            }

            outputModifications(pw, rs.getInt("id"), indent);         
        }
        pw.print(indent + "<MS_Data> <!-- looping element -->\n");
        outputDataFiles(pw, experimentId, indent); 
        pw.print(indent + "</MS_Data>\n");
        return outputProteinHits(pw, experimentId, indent);
      
    }
    // also output the parameterFile
    private boolean getStaticMod(StringBuffer sb, int searchParamId, String indent) throws SQLException, IOException {
        boolean hasResult = false;
        String q = "SELECT parameters from db_search_parameters where id = " + searchParamId;
        ResultSet rs = executeQuery(q);
        String modBTag = indent + "<Modification>\n";
        String modETag = indent + "</Modification>\n";
        indent += "\t";
        if(rs.next()) {
            String params = rs.getString("parameters");
            BufferedReader br = new BufferedReader(new StringReader(params)); 
            String line = br.readLine();
            while (line != null) {
                if(line.startsWith("add_")) {
                    String [] splitBySemiColon = line.split(";");
                    String [] splitByEqualSign = splitBySemiColon[0].trim().split("=");
                    double massShift = Double.parseDouble(splitByEqualSign[1].trim());
                    if(massShift != 0) {
                        hasResult = true;
                        //System.out.println(splitByEqualSign[0] + splitByEqualSign[1]);
                        sb.append(modBTag);
                        String residue = splitByEqualSign[0].split("_")[1];

                        //System.out.println("static mod residue: " + residue + "\tmassShift: " + massShift);
                        sb.append(outputContent(indent, "Modified_Residue", residue));
                        sb.append(outputContent(indent, "Mass_Shift", ""+massShift));
                        sb.append(outputContent(indent, "Is_Diff_Mod", "false"));

                        sb.append(modETag);
                    }
                }

                line = br.readLine();
            }
            br.close();
        }
        return hasResult;
    }
    private void outputModifications(PrintWriter pw, int searchParamId, String indent) throws SQLException, IOException {
        StringBuffer sb = new StringBuffer();
        String q = "SELECT mt.* FROM modification_type mt, db_search_parameters_has_modification_type ds " +
                   "WHERE ds.db_search_parameters_id = " + searchParamId +  
                   " AND mt.id = modification_type_id";
        ResultSet rs = executeQuery(q);

        boolean hasResult = false;
        sb.append(indent + "<Modifications>\n");
        String subIndent = indent + "\t";
        String moreIndent = "\t" + subIndent;
 
        hasResult = getStaticMod(sb, searchParamId, subIndent);

        while(rs.next()) {
            hasResult = true;
            sb.append(subIndent + "<Modification>\n");
            sb.append(outputContent(moreIndent, "Modified_Residue", ""+rs.getString("modified_residue")));
            sb.append(outputContent(moreIndent, "Mass_Shift", ""+rs.getFloat("mass_shift")));
            sb.append(outputContent(moreIndent, "Is_Diff_Mod", "true"));
            sb.append(subIndent + "</Modification>\n");
        }
        sb.append(indent + "</Modifications>\n");
        if(hasResult) {
            pw.print(sb.toString());
        }
    }
    private int outputSubExperiment(PrintWriter pw, int experimentId, String indent) throws SQLException, IOException {
        indent += "\t";
        pw.print(outputContent(indent, "Mass_Spec_Machine_ID", experimentId2InstrumentType(experimentId)));        
        pw.print(indent + "<DB_Search> <!-- looping element -->\n");    
        int numProteins = outputDbsearch(pw, (experimentId), indent);    
        pw.print(indent + "</DB_Search>\n");    
        return numProteins;
    }
    // temporary solution
    private String experimentId2GenderId(int exp_id) throws SQLException {

        String query = "SELECT donor_id FROM sample s, msp_experiment m WHERE project_id = 16 and sample_id = s.id and m.id = " + exp_id; 
        ResultSet rs = executeQuery(query);
        while(rs.next()) {
            String donorId = rs.getString("donor_id");
            if (donorId == null ) {
                return "O";
            } else if(donorId.startsWith("3")) {
                return "F";
            } else {
                return "M";
            }
        }
        return "O";
    }
    private int outputExperiment(PrintWriter pw, int experimentId, String indent) throws SQLException, IOException {
        pw.print("\t<Experiment> <!-- looping element -->\n");
        indent += "\t";
        pw.print(outputContent(indent, "Experiment_Submission_ID", "HSPP_JY_" + experimentId));
        pw.print(outputContent(indent, "Laboratory_ID", "1"));
        pw.print(outputContent(indent, "Original_Experiment_ID", experimentId + ""));
        String sampleSourceId = experimentId2SampleSourceId(experimentId);
        pw.print(outputContent(indent, "Sample_Source_ID", sampleSourceId));
        // need to add Gender_ID here
        pw.print(outputContent(indent, "Gender_ID", experimentId2GenderId(experimentId)));

        pw.print(indent + "<Sub_Experiment> <!-- looping element -->\n");
        int numProteins = outputSubExperiment(pw, experimentId, indent);
        pw.print(indent + "</Sub_Experiment>\n");
        pw.print("\t</Experiment>\n"); 
        return numProteins;
    }
    public void experimentIds2Xml(int [] ids) throws SQLException, IOException {
        
        xmlFile.print(HEADER);
        for(int id : ids) {
            if(!experimentHasResults(id)) {
                continue;
            }
            String indent = "\t";
            int numProteins = outputExperiment(xmlFile, id, indent); 
            printLog("Experiment_Submission_ID: HSPP_JY_" + id + "\texperiment_id: " + id + "\tnumProteins: " + numProteins);
             
        }
        xmlFile.print(TAIL);
        xmlFile.close();
    } 
    public static void main(String [] args) {
        try {
            Pmsdb pdb = new Pmsdb("newpmsdb", args[0]);
            Date d = new Date(System.currentTimeMillis());
            pdb.printLog(d.toString());
            pdb.printLog("Experiments outputed: ");   
          /* 
            int [] projectIds = pdb.projectId2ExperimentIds(16, 434);
            for(int id : projectIds) {
                System.out.println("project id: " + id);
            }
          */
          /*  
            //for(int id : pdb.projectId2ExperimentIds(16)) {
            //for(int id : pdb.projectId2ExperimentIds(16, 420)) {
            for(int id : pdb.projectId2ExperimentIds(16, 678)) {
                System.out.println(id + "\t" + pdb.experimentHasResults(id));
            }

            System.out.println("exiting");
           */
            //int [] experimentIds = pdb.projectId2ExperimentIds(16, 420);
            //int [] experimentIds = {139, 140};
            //int [] experimentIds = {892, 894}; // for two-pass results
            //int [] experimentIds = {827};
            int [] experimentIds = {892};
            pdb.experimentIds2Xml(experimentIds);
            pdb.printLog("Total number of unique proteins submited: " + pdb.proteins.size());
            pdb.printLog("XML file generated: " + pdb.xmlFileName + "\n\n\n");
            
            pdb.closeLog();
            pdb.closeConnection();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }

    /**
     * @fn void closeConnection()
     * Close the DB connection
     * @return void
     * @exception SQLException
     */
    public void closeConnection() throws SQLException {
        super.closeConnection();
    }

}






