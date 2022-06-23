package scripts.imscore;

import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

public class ShuffleFasta {

    public static void main(String args[]) throws IOException {
        //String fname = "/home/ip2/data/test.fasta";
        String fname = args[0];
        for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fname)); itr.hasNext(); ) {
            Fasta fasta = itr.next();
            String defline = fasta.getDefline();

            System.out.println(">Reverse_" + defline);
            //System.out.println(fasta.getSequence());
            System.out.println(shuffleString(fasta.getSequence()));
        }
    }

    public static String shuffleString(String string) {
        List<String> letters = Arrays.asList(string.split(""));
        Collections.shuffle(letters);
        StringBuilder builder = new StringBuilder();
        for (String letter : letters) {
            builder.append(letter);
        }
        return builder.toString();
    }
}
