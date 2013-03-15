package org.broadinstitute;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broadinstitute.pairhmm.*;
import org.broadinstitute.utils.QualityUtils;

import java.io.*;
import java.util.*;


/**
 * Evaluating routine for the PairHMM
 *
 * @author Mauricio Carneiro
 * @since 8/15/12
 */
public class Evaluate {

    public final List<PairHMM> pairHMM = new ArrayList<PairHMM>(10);
    public static boolean useRootbeer;

    private static final int X_METRIC_LENGTH = 10000;
    private static final int Y_METRIC_LENGTH = 10000;

    private static Logger logger = Logger.getLogger("Main");

    public Evaluate(Set<String> args) {
        useRootbeer = false;
        final boolean debug = args.contains("-d") || args.contains("--debug");
        logger.setLevel(debug ? Level.DEBUG : Level.INFO);
        BasicConfigurator.configure();
    }

    public void initializeHMMs(Set<String> args) {
        boolean addedHMMs = false;

        if (args.contains("--all") || args.contains("--approximate")) {
            logger.info("Initializing ApproximatePairHMM");
            pairHMM.add(new ApproximatePairHMM());
            addedHMMs = true;
        }
        if (args.contains("--all") || args.contains("--standard")) {
            logger.info("Initializing StandardPairHMM");
            pairHMM.add(new StandardPairHMM());
            addedHMMs = true;
        }
        if (args.contains("--all") || args.contains("--exact")) {
            logger.info("Initializing ExactPairHMM");
            pairHMM.add(new ExactPairHMM());
            addedHMMs = true;
        }
        if (args.contains("--all") || args.contains("--max")) {
            logger.info("Initializing MaxPairHMM");
            pairHMM.add(new MaxPairHMM());
            addedHMMs = true;
        }
        if (args.contains("--all") || args.contains("--experimental")) {
            logger.info("Initializing ExperimentalPairHMM");
            pairHMM.add(new ExperimentalPairHMM());
            addedHMMs = true;
        }
        if (!addedHMMs || args.contains("--all") || args.contains("--logless")) {
            logger.info("Initializing LoglessPairHMM");
            pairHMM.add(new LoglessPairHMM());
        }
        if(args.contains("--all") || args.contains("--rootbeer")){
            useRootbeer = true;
        }
    }

    /**
     * Initializes and computes the Pair HMM matrix.
     * <p/>
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
    public double runhmm(PairHMM hmm, final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, final int hapStartIndex, final boolean recacheReadValues) {
        return hmm.subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, cleanupQualityScores(readQuals), insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);
    }

    private static String createFileName(String hmmName) throws IOException {
        int runNumber = 0;
        String filename;
        File testFileName;
        do {
            runNumber++;
            filename = hmmName + "." + runNumber + ".out";
            testFileName = new File(filename);
        } while (testFileName.exists());
        return filename;
    }

    /**
     * Ensures that all the qual scores have valid values
     *
     * @param readQuals read qualities array (modified in place)
     * @return the readQuals array for convenience
     */
    private byte[] cleanupQualityScores(byte[] readQuals) {
        for (int i = 0; i < readQuals.length; i++) {
            readQuals[i] = (readQuals[i] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[i] > Byte.MAX_VALUE ? Byte.MAX_VALUE : readQuals[i]));
        }
        return readQuals;
    }

    private void runTests(final PairHMM hmm, final Iterator<TestRow> testCache, String hmmName, String testSet) throws IOException {
        long totalTime = 0L;
        final String runName = hmmName + "." + testSet;
        final String filename = createFileName(runName);
        final FileWriter out = new FileWriter(filename);
        hmm.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);
        while (testCache.hasNext()) {
            final TestRow currentTest = testCache.next();
            final long startTime = System.nanoTime();
            final double likelihood = runhmm(hmm, currentTest.haplotypeBases, currentTest.readBases, currentTest.readQuals, currentTest.readInsQuals, currentTest.readDelQuals, currentTest.overallGCP, currentTest.haplotypeStart, currentTest.reachedReadValue);
            totalTime += System.nanoTime() - startTime;
            out.write("" + likelihood + "\n");
        }
        logger.info(String.format("%s test completed in %.3f secs, results written to: %s", hmmName, totalTime/1000000000.0, filename));
        out.close();
    }


    private static Iterator<TestRow> createIteratorFor(String filePath) throws FileNotFoundException {
        return new ParserIterator(filePath);
    }

    public static void main(final String[] argv) throws IOException {

        Set<String> args = new HashSet<String>(Arrays.asList(argv));
        if (args.size() < 1) {
            throw new RuntimeException("\r\nYou must specify a file name for input.\n" + "filename \n ----------------------------\n" + "Run with -d or --debug for debug output");
        } else {
            final Evaluate evaluate = new Evaluate(args);
            final List<String> testFiles = new LinkedList<String>();

            for (String arg : args) {
                if (!arg.startsWith("-")) {
                    logger.info("Adding test dataset: " + arg);
                    testFiles.add(arg);
                }
            }

            for (final String testSet : testFiles) {
                logger.info("Using " + testSet + " tests");

                final String[] testSplit = testSet.split("/");                           // get rid of the file path (if any)
                final String testName = testSplit[testSplit.length - 1].split("\\.")[0]; // get rid of the file extension
                    
                evaluate.initializeHMMs(args);
                final Iterator<PairHMM> pairHMMIterator = evaluate.pairHMM.iterator();
                while(pairHMMIterator.hasNext()) {
                    final PairHMM hmm = pairHMMIterator.next();
                    final String hmmName = hmm.getClass().getSimpleName();
                    logger.info("Running " + hmmName);
                    evaluate.runTests(hmm, createIteratorFor(testSet), hmmName, testName);
                    pairHMMIterator.remove();
                }

                if(useRootbeer){
                    logger.info("Running Rootbeer");
                    RootbeerEvaluate rb_evaluator = new RootbeerEvaluate();
                    final String hmmName = "RootbeerHmm";
                    final String runName = hmmName + "." + testSet;
                    rb_evaluator.runTests(createIteratorFor(testSet), hmmName, testName, createFileName(runName));
                }

                logger.info("Finished all HMMs for " + testSet + " tests");
            }
        }

    }

    /**
     * Classes and Functions used to parse and apply the tests to the
     * PairHMM class
     */

    static class ParserIterator implements Iterator<TestRow> {
        private XReadLines m_reader;
        private TestRow m_nextEntry;

        public ParserIterator(String filePath) throws FileNotFoundException {
            m_reader = new XReadLines(new File(filePath));
        }

        @Override
        public boolean hasNext() {
            if (!m_reader.hasNext()) {
                return false;
            }

            String line = m_reader.next();
            String[] p = line.split(" ");
            final boolean reachedRead = p[7].trim().equals("true");
            m_nextEntry = new TestRow(p[0].getBytes(), p[1].getBytes(), convertToByteArray(p[2]), convertToByteArray(p[3]), convertToByteArray(p[4]), convertToByteArray(p[5]), Integer.parseInt(p[6]), reachedRead);

            logger.debug("Haplotype Bases:" + p[0]);
            logger.debug("Read Bases:" + p[1]);
            logger.debug("Read Quals: " + p[2]);
            logger.debug("Read Ins Quals: " + p[3]);
            logger.debug("Read Del Quals: " + p[4]);
            logger.debug("Overall GCP:" + p[5]);
            logger.debug("Haplotype start:" + p[6]);
            logger.debug("Boolean (reached read values):" + p[7]);

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
        byte[] haplotypeBases;
        byte[] readBases;
        byte[] readQuals;
        byte[] readInsQuals;
        byte[] readDelQuals;
        byte[] overallGCP;
        int haplotypeStart;
        boolean reachedReadValue;

        TestRow(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] readInsQuals, byte[] readDelQuals, byte[] overallGCP, int haplotypeStart, boolean reachedReadValue) {
            this.haplotypeBases = haplotypeBases;
            this.readBases = readBases;
            this.readQuals = readQuals;
            this.readInsQuals = readInsQuals;
            this.readDelQuals = readDelQuals;
            this.overallGCP = overallGCP;
            this.haplotypeStart = haplotypeStart;
            this.reachedReadValue = reachedReadValue;
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
