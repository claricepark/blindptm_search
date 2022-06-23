package scripts.imscore;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class DTASelectTxtParser {

    public static void main(String args[]) throws IOException
    {

        run(args[0]);
        //System.out.println("aaa");
    }

    /**********************
     * Count number of PSM after given FDR
     */




    public static void run(String filename) throws IOException {

//        System.out.println(Double.parseDouble("9.0E-4"));
//        if(true) return;

         BufferedReader br = new BufferedReader(new FileReader(filename));
        String eachline;

//        int count=0;
        boolean isTarget=true;

        List<Peptide> list = new ArrayList<>();
        while (null != (eachline = br.readLine())) {
            //if(count++>50)
              //  break;
           // System.out.println(eachline);

            if(eachline.startsWith("H\t"))
                continue;

            if(!eachline.startsWith("D\t") && !eachline.startsWith("L\t") )
                continue;

            if(eachline.startsWith("L\t")) {
                String protein = eachline;
                if(protein.startsWith("L\tReverse"))
                    isTarget = false;
                else
                    isTarget = true;

                continue;
            }

            String[] arr = eachline.split("\t");
            Peptide pep = new Peptide(Double.parseDouble(arr[5]), isTarget);

            System.out.println(isTarget +"\t" + arr[5] + "\t" + eachline);

            list.add(pep);

          //  System.out.println(arr[5] + "\t" + eachline);
        }
//
////        System.out.println(list);
//        for(int i=0;i<20;i++)
//            System.out.println(list.get(i).getScore() + " " + list.get(i).isTarget());

        Collections.sort(list, new PeptideCompare());

//        System.out.println("==================");


     //   System.out.println(list);
        int targetCount=0;
        int falseCount=0;

        //for(int i=0;i<20;i++) {
        for(int i=0;i<list.size();i++) {
            double score = list.get(i).getScore();
            boolean target = list.get(i).isTarget();
//            System.out.println( + " " + );
            if(target)
                targetCount++;
            else
                falseCount++;

            double fdr =  falseCount/(double)(targetCount+falseCount);

            if(fdr>0.01)
                break;

        }

        System.out.println(targetCount + "\t" + falseCount + "\t" + (targetCount + falseCount) + "\t" + falseCount/(double)(targetCount+falseCount));

        /*
        ArrayList<Student> ar = new ArrayList<Student>();
        ar.add(new Student(111, "bbbb", "london"));
        ar.add(new Student(131, "aaaa", "nyc"));
        ar.add(new Student(121, "cccc", "jaipur"));

        System.out.println("Unsorted");
        for (int i = 0; i < ar.size(); i++)
            System.out.println(ar.get(i));

        Collections.sort(ar, new Sortbyroll());

        System.out.println("\nSorted by rollno");
        for (int i = 0; i < ar.size(); i++)
            System.out.println(ar.get(i));

        Collections.sort(ar, new Sortbyname());

        System.out.println("\nSorted by name");
        for (int i = 0; i < ar.size(); i++)
            System.out.println(ar.get(i));
    }
}
         */
        br.close();
    }

}

class Peptide {
    double score;
    boolean target;

    // Constructor
    public Peptide(double score, boolean target)
    {
        this.score = score;
        this.target = target;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public boolean isTarget() {
        return target;
    }

    public void setTarget(boolean target) {
        this.target = target;
    }
}

class PeptideCompare implements Comparator<Peptide> {
    // Used for sorting in ascending order of
    // roll number
    public int compare(Peptide p1, Peptide p2)
    {

        double diff = p2.score - p1.score;

        if(diff>0)
            return 1;
        else if (diff<0)
            return -1;
        else
            return 0;
//        return (int)(diff*10000);
    }
}
