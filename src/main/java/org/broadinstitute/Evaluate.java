package org.broadinstitute;

import org.broadinstitute.pairhmm.LoglessCachingPairHMM;
import org.broadinstitute.pairhmm.PairHMM;
import org.broadinstitute.utils.QualityUtils;

import java.io.*;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
  * Evaluating routine for the PairHMM
  *
  * @author Mauricio Carneiro
  * @since 8/15/12
  */
public class Evaluate {

    private static final long startTime = System.currentTimeMillis();
    private static final PairHMM pairHMM = new LoglessCachingPairHMM();
    private static final double PRECISION = 0.001;

    private static final int X_METRIC_LENGTH = 10000;
    private static final int Y_METRIC_LENGTH = 10000;

    public Evaluate() {
    }

    /**
      * Initializes and computes the Pair HMM matrix.
      *
      * Use this method if you're calculating the entire matrix from scratch.
      *
      * @param haplotypeBases reference sequence bases
      * @param readBases      comparison haplotype bases
      * @param readQuals      comparison haplotype base quals (phred-scaled)
      * @param insertionGOP   comparison haplotype insertion quals (phred-scaled)
      * @param deletionGOP    comparison haplotype deletion quals (phred-scaled)
      * @param overallGCP     comparison haplotype gap continuation quals (phred-scaled)
      * @return the likelihood of the alignment between read and haplotype
      */
    public static double hmm(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final byte[] insertionGOP,
                              final byte[] deletionGOP, final byte[] overallGCP,
                              final int hapStartIndex, final boolean recacheReadValues) {
        // ensure that all the qual scores have valid values
        for (int i = 0; i < readQuals.length; i++) {
            readQuals[i] = (readQuals[i] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[i] > Byte.MAX_VALUE ? Byte.MAX_VALUE : readQuals[i]));
        }

        return pairHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);
    }

    private static boolean runTests(Iterator<TestRow> testCache, boolean debug, boolean write) {

        pairHMM.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);

        if (write) {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output.txt"), false));
                while (testCache.hasNext()) {
                    TestRow currentTest = testCache.next();

                    Double result = hmm(currentTest.getHaplotypeBases(), currentTest.getReadBases(),
                            currentTest.getReadQuals(), currentTest.getReadInsQuals(),
                            currentTest.getReadDelQuals(), currentTest.getOverallGCP(),
                            currentTest.getHaplotypeStart(), currentTest.getReachedReadValye());

                    if (debug) {
                        System.out.printf(" Result:%4.3f%n" +
                                "==========================================================%n", result);
                    }
                    if (Math.abs(currentTest.getLikelihood() - result) > PRECISION) {
                        System.out.println("Wrong result. Expected " + currentTest.getLikelihood() + " , actual: " + result);
                        return false;
                    }
                    bw.write(String.format("%f", result));
                    bw.newLine();
                }
                bw.close();
                System.out.printf("%d - First run-through complete.%n", System.currentTimeMillis() - startTime);

            }  catch (Exception e) {//Catch exception if any
                throw new RuntimeException(e);
                }
        }    else   {
            while (testCache.hasNext()) {
                TestRow currentTest = testCache.next();

                Double result = hmm(currentTest.getHaplotypeBases(), currentTest.getReadBases(),
                        currentTest.getReadQuals(), currentTest.getReadInsQuals(),
                        currentTest.getReadDelQuals(), currentTest.getOverallGCP(),
                        currentTest.getHaplotypeStart(), currentTest.getReachedReadValye());

                if (debug) {
                    System.out.printf(" Result:%4.3f%n" +
                            "==========================================================%n", result);
                }
                if ((currentTest.getLikelihood() - result) > PRECISION) {
                    System.out.println("Wrong result. Expected " + currentTest.getLikelihood() + " , actual: " + result);
                    return false;
                }
            }
            System.out.printf("%d - No I/O run-through complete.%n", System.currentTimeMillis() - startTime);
        }
        return true;
    }


    private static Iterator<TestRow> parseFile(String filePath, boolean debug) {
        try {
            ParserIterator ret = new ParserIterator(filePath, debug);
            return ret;
        } catch (FileNotFoundException ex) {
            ex.printStackTrace(System.out);
            System.exit(0);
            return null;
        }
    }

    public static void main(String[] argv) {
        boolean debug = false;


        List<String> args = new LinkedList<String>(Arrays.asList(argv));
        if (args.contains("-debug")) {
            debug = true;
            args.remove("-debug");
        }
        if (args.size() < 1) {
            throw new RuntimeException("\r\nYou must specify a file name for input.Format below:\n" +
                    "filename \n ----------------------------\n" +
                    "Run with -debug for debug output");
        } else {
            for (String arg : args) {
                System.out.println("testing file: " + arg);
                Iterator<TestRow> iter = parseFile(arg, debug);
                if (runTests(iter, debug, true)) {
                    System.out.println("Tests successful");
                } else {
                    System.out.println("Tests unsuccessful");
                }
            }
        }

    }

    /*************************************************************************
     * Classes and Functions used to parse and apply the tests to the
     * PairHMM class
     ************************************************************************/

    static class ParserIterator implements Iterator<TestRow> {

      private XReadLines m_reader;
      private TestRow m_nextEntry;
      private boolean m_debug;

      public ParserIterator(String filePath, boolean debug) throws FileNotFoundException {
        m_reader = new XReadLines(new File(filePath));
        m_debug = debug;
      }

      @Override
      public boolean hasNext() {
        if(m_reader.hasNext() == false){
          return false;
        }

        String line = m_reader.next();
        String[] p = line.split(" ");
        boolean reachedRead;
        if (p[7].trim().equals("true")) {
          reachedRead = true;
        } else {
          reachedRead = false;
        }

        TestRow testRow = new TestRow(p[0].getBytes(), p[1].getBytes(),
          convertToByteArray(p[2]), convertToByteArray(p[3]),
          convertToByteArray(p[4]), convertToByteArray(p[5]),
          Integer.parseInt(p[6]), reachedRead,
          Double.parseDouble(p[8]));

        m_nextEntry = testRow;

        if (m_debug) {
          System.out.println("Haplotype Bases:" + p[0]);
          System.out.println("Read Bases:" + p[1]);
          System.out.println("Read Quals: " + p[2]);
          System.out.println("Read Ins Quals: " + p[3]);
          System.out.println("Read Del Quals: " + p[4]);
          System.out.println("Overall GCP:" + p[5]);
          System.out.println("Haplotype start:" + p[6]);
          System.out.println("Boolean (reached read values):" + p[7]);
          System.out.println("result, likelihood:" + p[8]);
        }
        return true;
      }

      private byte[] convertToByteArray(String quals) {
        byte[] output = quals.getBytes();

        for (int i = 0; i < output.length; i++) {
          output[i] -= 33;
        }

        return output;
      }

      @Override
      public TestRow next() {
        return m_nextEntry;
      }

      @Override
      public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
      }

    }

    static class TestRow {
        private byte[] haplotypeBases;
        private byte[] readBases;
        private byte[] readQuals;
        private byte[] readInsQuals;
        private byte[] readDelQuals;
        private byte[] overallGCP;
        private int haplotypeStart;
        private boolean reachedReadValue;
        private Double likelihood;

        TestRow(byte[] haplotypeBases, byte[] readBases,
                byte[] readQuals, byte[] readInsQuals,
                byte[] readDelQuals, byte[] overallGCP,
                int haplotypeStart, boolean reachedReadValue,
                Double likelihood) {

            this.haplotypeBases = haplotypeBases;
            this.readBases = readBases;
            this.readQuals = readQuals;
            this.readInsQuals = readInsQuals;
            this.readDelQuals = readDelQuals;
            this.overallGCP = overallGCP;
            this.haplotypeStart = haplotypeStart;
            this.reachedReadValue = reachedReadValue;
            this.likelihood = likelihood;
        }

        public byte[] getHaplotypeBases() {
            return this.haplotypeBases;
        }

        public byte[] getReadBases() {
            return this.readBases;
        }

        public byte[] getReadQuals() {
            return this.readQuals;
        }

        public byte[] getReadInsQuals() {
            return this.readInsQuals;
        }

        public byte[] getReadDelQuals() {
            return this.readDelQuals;
        }

        public byte[] getOverallGCP() {
            return this.overallGCP;
        }

        public int getHaplotypeStart() {
            return this.haplotypeStart;
        }

        public boolean getReachedReadValye() {
            return this.reachedReadValue;
        }

        public double getLikelihood() {
            return this.likelihood;
        }

    }

    static class XReadLines implements Iterator<String>, Iterable<String> {
        private final BufferedReader in;      // The stream we're reading from
        private String nextLine = null;       // Return value of next call to next()
        private final boolean trimWhitespace;
        private final String commentPrefix;

        public XReadLines(final File filename) throws FileNotFoundException {
            this(new FileReader(filename), true, null);
        }

        /**
         * Creates a new xReadLines object to read lines from an bufferedReader
         *
         * @param reader         file name
         * @param trimWhitespace trim whitespace
         * @param commentPrefix  prefix for comments or null if no prefix is set
         */
        public XReadLines(final Reader reader, final boolean trimWhitespace, final String commentPrefix) {
            this.in = (reader instanceof BufferedReader) ? (BufferedReader) reader : new BufferedReader(reader);
            this.trimWhitespace = trimWhitespace;
            this.commentPrefix = commentPrefix;
            try {
                this.nextLine = readNextLine();
            } catch (IOException e) {
                throw new IllegalArgumentException(e);
            }
        }

        /**
         * I'm an iterator too...
         *
         * @return an iterator
         */
        public Iterator<String> iterator() {
            return this;
        }

        public boolean hasNext() {
            return this.nextLine != null;
        }

        /**
         * Actually reads the next line from the stream, not accessible publicly
         *
         * @return the next line or null
         * @throws java.io.IOException if an error occurs
         */
        private String readNextLine() throws IOException {
            String nextLine;
            while ((nextLine = this.in.readLine()) != null) {
                if (this.trimWhitespace) {
                    nextLine = nextLine.trim();
                    if (nextLine.length() == 0)
                        continue;
                }
                if (this.commentPrefix != null)
                    if (nextLine.startsWith(this.commentPrefix))
                        continue;
                break;
            }
            return nextLine;
        }

        /**
         * Returns the next line (optionally minus whitespace)
         *
         * @return the next line
         */
        public String next() {
            try {
                String result = this.nextLine;
                this.nextLine = readNextLine();

                // If we haven't reached EOF yet
                if (this.nextLine == null) {
                    in.close();             // And close on EOF
                }

                // Return the line we read last time through.
                return result;
            } catch (IOException e) {
                throw new IllegalArgumentException(e);
            }
        }

        // The file is read-only; we don't allow lines to be removed.
        public void remove() {
            throw new UnsupportedOperationException();
        }

    }
}
