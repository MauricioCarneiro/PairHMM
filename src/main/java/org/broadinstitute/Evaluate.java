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

    private static final double PRECISION = 0.001;
    private static final int X_METRIC_LENGTH = 10000;
    private static final int Y_METRIC_LENGTH = 10000;

    private static Logger logger = Logger.getLogger("Main");

    public Evaluate(Set<String> args) {
        boolean addedHMMs = false;
        useRootbeer = false;

        // Setup the logger
        final boolean debug = args.contains("-d") || args.contains("--debug");
        logger.setLevel(debug ? Level.DEBUG : Level.INFO);
        BasicConfigurator.configure();

        if (args.contains("--all") || args.contains("--caching")) {
            logger.info("Including CachingPairHMM");
            pairHMM.add(new CachingPairHMM());
            addedHMMs = true;
        }
        if (args.contains("--all") || args.contains("--original")) {
            logger.info("Including OriginalPairHMM");
            pairHMM.add(new OriginalPairHMM());
            addedHMMs = true;
        }
        if (args.contains("--all") || args.contains("--exact")) {
            logger.info("Including ExactPairHMM");
            pairHMM.add(new ExactPairHMM());
            addedHMMs = true;
        }
        if (!addedHMMs || args.contains("--all") || args.contains("--logless")) {
            logger.info("Including LoglessCachingPairHMM");
            pairHMM.add(new LoglessCachingPairHMM());
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

    private void runTests(final PairHMM hmm, final Iterator<TestRow> testCache) {
        final long startTime = System.currentTimeMillis();

        boolean testsPassed = true;
        hmm.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);

        while (testCache.hasNext()) {
            TestRow currentTest = testCache.next();

            Double result = runhmm(hmm, currentTest.getHaplotypeBases(), currentTest.getReadBases(), currentTest.getReadQuals(), currentTest.getReadInsQuals(), currentTest.getReadDelQuals(), currentTest.getOverallGCP(), currentTest.getHaplotypeStart(), currentTest.getReachedReadValye());

            logger.debug(String.format(" Result:%4.3f",result));
            logger.debug("==========================================================%n");
            if (Math.abs(currentTest.getLikelihood() - result) > PRECISION) {
                logger.error("Wrong result. Expected " + currentTest.getLikelihood() + " , actual: " + result);
                testsPassed = false;
            }
        }
        if (testsPassed)
            logger.info(String.format("All tests PASSED in %.3f secs", (double) (System.currentTimeMillis() - startTime)/1000));

        hmm.clear();
    }


    private static Iterator<TestRow> createIteratorFor(String filePath) throws FileNotFoundException {
        return new ParserIterator(filePath);
    }

    public static void main(final String[] argv) throws FileNotFoundException {
        Set<String> args = new HashSet<String>(Arrays.asList(argv));
        if (args.size() < 1) {
            throw new RuntimeException("\r\nYou must specify a file name for input.\n" + "filename \n ----------------------------\n" + "Run with -d or --debug for debug output");
        } else {
            final Evaluate evaluate = new Evaluate(args);
            final List<String> testFiles = new LinkedList<String>();

            for (String arg : args) {
                if (!arg.startsWith("-")) {
                    logger.info("adding testset: " + arg);
                    testFiles.add(arg);
                }
            }

            for (final String testSet : testFiles) {
                logger.info("RUNNING " + testSet + " tests");
                for(final PairHMM hmm : evaluate.pairHMM) {
                    logger.info("EVALUATING " + hmm.getClass().getSimpleName());
                    evaluate.runTests(hmm, createIteratorFor(testSet)); // runTests will output the detailed test results
                }
                if(useRootbeer){
                    logger.info("EVALUATING Rootbeer");
                    RootbeerEvaluate rb_evaluator = new RootbeerEvaluate(PRECISION);
                    rb_evaluator.runTests(createIteratorFor(testSet));
                }
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
            m_nextEntry = new TestRow(p[0].getBytes(), p[1].getBytes(), convertToByteArray(p[2]), convertToByteArray(p[3]), convertToByteArray(p[4]), convertToByteArray(p[5]), Integer.parseInt(p[6]), reachedRead, Double.parseDouble(p[8]));

            logger.debug("Haplotype Bases:" + p[0]);
            logger.debug("Read Bases:" + p[1]);
            logger.debug("Read Quals: " + p[2]);
            logger.debug("Read Ins Quals: " + p[3]);
            logger.debug("Read Del Quals: " + p[4]);
            logger.debug("Overall GCP:" + p[5]);
            logger.debug("Haplotype start:" + p[6]);
            logger.debug("Boolean (reached read values):" + p[7]);
            logger.debug("result, likelihood:" + p[8]);

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

    public static class TestRow {
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
