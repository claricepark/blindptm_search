/**
 * @file FileUtils.java
 * This is the source file for edu.scripps.pms.util.spectrum.FileUtils
 * @author Tao Xu
 * @date $Date
 */



package edu.scripps.pms.util.io;

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;

public class FileUtils {
    public static ArrayList<String> getFilePrefixes(String dir, String suffix) throws IOException {

         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(suffix)) {
                 //System.out.println("got ms2 file: " + s);
                 files.add(s.split(suffix)[0]);
             }
         }
         return files;

    }
    // get all files in pwd that ends with suffix
    public static ArrayList<String> getFiles(String dir, String suffix) throws IOException {

         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(suffix)) {
                 //System.out.println("got ms2 file: " + s);
                 files.add(s);
             }
         }
         return files;

    }
    // get all files in pwd that ends with one of the suffixes
    public static ArrayList<String> getFiles(String dir, ArrayList<String> suffixes) throws IOException {

         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             for(String suffix : suffixes) {
                if(s.endsWith(suffix)) {
                    //System.out.println("got ms2 file: " + s);
                    files.add(s);
                }
             }
         }
         return files;

    }
}
