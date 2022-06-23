package clariceparky.informatics;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;

public class ReadWriteFile {
    public static void main(String[] args) throws Exception {
//google java type: integer, double, String, etc.
        String fileName = "/Users/cpark/Downloads/ms2.txt";
        String output = "/Users/cpark/Downloads/ms2_out.txt";

        BufferedOutputStream bo =
                new BufferedOutputStream(new FileOutputStream(output));

        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = br.readLine();
        while(line != null) {
            System.out.println(line);
            line = br.readLine();
            //br.write(string);
        }

        br.close();
        bo.close();


    }
}
