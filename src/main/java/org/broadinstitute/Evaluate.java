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

    private static PairHMM pairHMM;
    private static final double PRECISION = 0.001;

    private static final int X_METRIC_LENGTH = 10000;
    private static final int Y_METRIC_LENGTH = 10000;

    private static Logger logger = Logger.getLogger("Main");

    public Evaluate(Set<String> args) {
        final boolean debug = args.contains("-d") || args.contains("--debug");
        logger.setLevel(debug ? Level.DEBUG : Level.INFO);
        BasicConfigurator.configure();
        if (args.contains("--caching")) {
            /*logger.info("Using CachingPairHMM");*/
            pairHMM = new CachingPairHMM();
        }
        else if (args.contains("--original")) {
            /*logger.info("Using OriginalPairHMM");*/
            pairHMM = new OriginalPairHMM();
        }
        else if (args.contains("--exact")) {
            /*logger.info("Using ExactPairHMM");*/
            pairHMM = new ExactPairHMM();
        }
        else {
            /*logger.info("Using LoglessCachingPairHMM");*/
            pairHMM = new LoglessCachingPairHMM();
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
    public double hmm(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, final int hapStartIndex, final boolean recacheReadValues) {
        return pairHMM.subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, cleanupQualityScores(readQuals), insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);
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


    private void runTestsImpa(Iterator<TestRow> testCache, String output_fn) {
        final long startTime = System.currentTimeMillis();
        long t1, comp_time = 0;

        pairHMM.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);

        try
        {
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_fn), false));
            while (testCache.hasNext())
            {
                TestRow currentTest = testCache.next();
                t1 = System.currentTimeMillis();
                Double result = hmm(currentTest.getHaplotypeBases(), currentTest.getReadBases(), currentTest.getReadQuals(), currentTest.getReadInsQuals(), currentTest.getReadDelQuals(), currentTest.getOverallGCP(), currentTest.getHaplotypeStart(), currentTest.getReachedReadValye());
                comp_time += (System.currentTimeMillis() - t1);
                bw.write(String.format("%f", result));
                bw.newLine();
            }
            bw.close();
            System.out.printf("COMPUTATION_TIME %d%n", comp_time);
        }  catch (Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    private void runTests(Iterator<TestRow> testCache) {
        final long startTime = System.currentTimeMillis();

        boolean testsPassed = true;
        pairHMM.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);

        while (testCache.hasNext()) {
            TestRow currentTest = testCache.next();

            Double result = hmm(currentTest.getHaplotypeBases(), currentTest.getReadBases(),
                    currentTest.getReadQuals(), currentTest.getReadInsQuals(),
                    currentTest.getReadDelQuals(), currentTest.getOverallGCP(),
                    currentTest.getHaplotypeStart(), currentTest.getReachedReadValye());

            logger.debug(String.format(" Result:%4.3f",result));
            logger.debug("==========================================================%n");
            if ((currentTest.getLikelihood() - result) > PRECISION) {
                logger.error("Wrong result. Expected " + currentTest.getLikelihood() + " , actual: " + result);
                testsPassed = false;
            }
        }
        if (testsPassed)
            logger.info(String.format("All tests PASSED in %.3f secs", (double) (System.currentTimeMillis() - startTime)/1000));
    }


    private static Iterator<TestRow> parseFile(String filePath) throws FileNotFoundException {
        return new ParserIterator(filePath);
    }

    public static void main(String[] argv) throws FileNotFoundException {
        Set<String> args = new HashSet<String>(Arrays.asList(argv));
        if (args.size() < 1) {
            throw new RuntimeException("\r\nYou must specify a file name for input.\n" + "filename \n ----------------------------\n" + "Run with -d or --debug for debug output");
        } else {
            Evaluate evaluate = new Evaluate(args);

			/*
				java -Xmx4g -jar build/libs/PairHMM-0.1.jar <input path> <output path> --impa-mode
			*/
			if (args.contains("--impa-mode")) {
				/*
                System.out.printf("%s%n", argv[0]);
                System.out.printf("%s%n", argv[1]);
                return;
				*/
				Iterator<TestRow> iter = parseFile(argv[0]);
				evaluate.runTestsImpa(iter, argv[1]);
			} else {
				for (String arg : args) {
					if (!arg.startsWith("-")) {
						logger.info("parsing file " + arg);
						Iterator<TestRow> iter = parseFile(arg);
						logger.info("running tests");
						evaluate.runTests(iter);
					}
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
