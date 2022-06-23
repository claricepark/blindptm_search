Thank you for licensing GutenTag for use in your lab.

Here's how to install it on your computer:


2) Extract the Zip into the directory C:\GutenTag on your hard disk.
   The password is "peptide" to open the zip.

3) Copy "GutenTag.bat" to your c:\windows or c:\winnt directory.

4) Test that your computer has Java installed: open a command prompt
   (usually Start -> All Programs -> Accessories -> MS-DOS command
   prompt), and type "java". If a help message shows up, Java is
   working on your computer. If you get an error message instead,
   you first need to install Java. You can get a JRE (Java Runtime
   Environment) at http://java.sun.com.


Here's how to use the program:

1) Create a directory on your hard disk for the test.

2) Copy the MS2 files you want to analyze to the new directory.

3) Copy the database you want to use to the directory (actually, the
   database just needs to be somewhere where Windows can see it).

4) Copy the files GutenTag.ini and Isotoper.ini to the test directory.

5) Modify the database path in GutenTag.ini to point to the database
   to be searched.

6) Open a command prompt (may be called "MS-DOS prompt" on some
   computers) and cd to the test directory.

7) Type "GutenTag" and hit enter.  Use "GutenTag --dta" instead if
   you're using DTA spectrum files instead of MS2 files.

GutenTag results are readable and reviewable by DTASelect.  In tryptic
peptides, a score of 12 or higher is likely to be a worthwhile match.
