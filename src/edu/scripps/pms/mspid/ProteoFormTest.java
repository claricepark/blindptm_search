package edu.scripps.pms.mspid;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSequenceWithPTMs;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.proteoform_dbindex.model.PTMCodeObj;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.masses.AssignMass;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Created by Titus Jung titusj@scripps.edu on 9/10/19.
 */
public class ProteoFormTest {

    public static void main(String [] args) throws DBIndexStoreException {
        // input parameters file
        File paramFile = new File("/home/yateslab/projectData/prolucid/3101ProteoformWork/0509SalvaNewFormat/blazmass_Q13523.params");
        String fasta = "/home/yateslab/projectData/prolucid/3101ProteoformWork/0509SalvaNewFormat/Q13523.fasta";
        int maxMisclevages = 3;
        double minPrcMass = 600;
        double maxPrcMass = 6000;
        String residues = "KR";
        String noCut = "P";
        int offset =0;
        boolean semiCleavage = false;



// maximum number of variations (sequence variations and PTMs) per peptide
        int maxNumVariationsPerPeptide = 4;

// Uniprot repository version release
// null for latest version or "2019_05" for May 2019 version, for example
        String uniprotVersion = null;
        File uniprotReleasesFolder = new File("/home/yateslab/projectData/prolucid/3101ProteoformWork/0509SalvaNewFormat/");

// Uniprot annotations retriever. It will retrieve the annotations to folder uniprotReleasesFolder
        UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
        //  System.out.println("TESTTTING ");
// Create ProteoformDBIndexInterface instance
        DBIndexSearchParams params =  DBIndexImpl.getDefaultDBIndexParamsForProteoformAnalysis(fasta,maxMisclevages,
                minPrcMass,maxPrcMass, residues,noCut,offset,false,null,"Reverse_",
                "/home/yateslab/projectData/prolucid/3101ProteoformWork/0509SalvaNewFormat/");
        ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(params, 4);

        //ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, uplr, uniprotVersion, maxNumVariationsPerPeptide);

        List<IndexedSequence> sequences = proteoformDBIndex.getSequences(2275.0, 20.0);
        System.out.println("PRINTING SEQ: ");

        ExtendedAssignMass mass = ExtendedAssignMass.getInstance(false,null);
        for(IndexedSequence indexedSequence: sequences)
        {
            if (indexedSequence instanceof IndexedSequenceWithPTMs) {
                final IndexedSequenceWithPTMs indexedSequenceWithPTM = (IndexedSequenceWithPTMs) indexedSequence;
             //   Pattern pattern = Pattern.compile("[A-Z][[+-0-9\.]*]");
               List<PTM> ptmSet =  indexedSequenceWithPTM.getPtms();
               for(PTM ptm : ptmSet)
               {
                   short code = ptm.getPtmCode();
                   PTMCodeObj obj =  mass.getPTMbyPTMCode( code);
                   System.out.println(">>><<<\t"+obj.getDescription()+"\t"+ptm.getPosInPeptideByte());
               }
               System.out.println(indexedSequenceWithPTM.getSequence() + "\t"
                        + IndexUtil.calculateMass(indexedSequenceWithPTM.getSequence()) + "\t"
                        + indexedSequenceWithPTM.getModSequence() + "\t" + indexedSequenceWithPTM.getMass());

            } else {
                System.out.println(indexedSequence.getSequence() + IndexUtil.calculateMass(indexedSequence.getSequence())
                        + "\t" + indexedSequence.getModSequence() + "\t" + +indexedSequence.getMass());
            }
        }
    }
}
