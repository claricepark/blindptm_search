package scripts.imscore;

import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.io.SQTParser;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.sqt.SQTPeptide;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;


public class DTASelectImScoreCompare {
    public static void main(String args[]) throws IOException
    {
        compareDTASelectImScore(args[0], args[1]);
    }

    public static void compareDTASelectImScore(String dtaselect, String sqtFile) throws IOException {

        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselect);

        Iterator<Protein> pitr = reader.getProteins();
        int count=0;

        Set<String> set = new HashSet<>();
        for (Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            Protein protein = itr.next();
            if(protein.getLocus().startsWith("Reve") || protein.getLocus().startsWith("contam"))
                continue;

            for(Iterator<Peptide> pepItr = protein.getPeptideList().iterator(); pepItr.hasNext(); )
            {
                Peptide peptide = pepItr.next();
               // System.out.println(peptide.getSequence());

                String[] arr = peptide.getPeptideLine();

                String key = peptide.getSequence() + peptide.getChargeState();
                if(set.contains(key))
                    continue;

                System.out.println(arr[12] + "\t" + arr[13] + "\t" + arr[14] + "\t" + peptide.getSequence() + "\t" + peptide.getChargeState());
                set.add(key);

               // peptide.getPredi
                count++;
            }
        }

        System.out.println(count);
    }

}
