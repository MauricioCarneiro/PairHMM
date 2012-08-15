/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */


import java.io.*;
import java.util.*;

/**
 * Hidden Markov Model for Haplotype assembly
 *
 * @author Mauricio Carneiro
 * @since 8/15/2012
 */

public class PairHMM {
    private static final int MAX_CACHED_QUAL = (int) Byte.MAX_VALUE;
    private static final byte DEFAULT_GOP = (byte) 45;
    private static final byte DEFAULT_GCP = (byte) 10;
    public final static byte MIN_USABLE_Q_SCORE = 6;
    private static double[] qualToErrorProbLog10Cache = new double[256];
    private static double[] qualToProbLog10Cache = new double[256];
    public static final double[] log10Cache;
    public static final double[] log10FactorialCache;
    private static final double[] jacobianLogTable;
    private static final double JACOBIAN_LOG_TABLE_STEP = 0.001;
    private static final double JACOBIAN_LOG_TABLE_INV_STEP = 1.0 / 0.001;
    private static final double MAX_JACOBIAN_TOLERANCE = 8.0;
    private static final int JACOBIAN_LOG_TABLE_SIZE = (int) (MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1;
    private static final int MAXN = 50000;
    private static final int LOG10_CACHE_SIZE = 4 * MAXN;  // we need to be able to go up to 2*(2N) when calculating some of the coefficients
    private static final long startTime = System.currentTimeMillis();

    static {
        for (int i = 0; i < 256; i++) qualToErrorProbLog10Cache[i] = qualToErrorProbLog10Raw(i);
        for (int i = 0; i < 256; i++) qualToProbLog10Cache[i] = qualToProbLog10Raw(i);

        log10Cache = new double[LOG10_CACHE_SIZE];
        log10FactorialCache = new double[LOG10_CACHE_SIZE];
        jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];

        log10Cache[0] = Double.NEGATIVE_INFINITY;
        for (int k = 1; k < LOG10_CACHE_SIZE; k++) {
            log10Cache[k] = Math.log10(k);
            log10FactorialCache[k] = log10FactorialCache[k - 1] + log10Cache[k];
        }

        for (int k = 0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
            jacobianLogTable[k] = Math.log10(1.0 + Math.pow(10.0, -((double) k) * JACOBIAN_LOG_TABLE_STEP));

        }
    }

    public static void initializeArrays(final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray,
                                        final int X_METRIC_LENGTH) {

        for (int iii = 0; iii < X_METRIC_LENGTH; iii++) {
            Arrays.fill(matchMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[iii], Double.NEGATIVE_INFINITY);
        }

        // the initial condition
        matchMetricArray[1][1] = 0.0; // Math.log10(1.0);

    }

    static private double qualToProbLog10Raw(int qual) {
        return Math.log10(1.0 - qualToErrorProbRaw(qual));
    }

    static public double qualToProbLog10(byte qual) {
        return qualToProbLog10Cache[(int) qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    /**
     * Convert a quality score to a probability of error.  This is the Phred-style
     * conversion, *not* the Illumina-style conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0 - 255)
     * @return a probability (0.0 - 1.0)
     */
    static private double qualToErrorProbRaw(int qual) {
        return qualToErrorProb((double) qual);
    }

    public static double qualToErrorProb(final double qual) {
        return Math.pow(10.0, ((double) qual) / -10.0);
    }

    static private double qualToErrorProbLog10Raw(int qual) {
        return ((double) qual) / -10.0;
    }

    static public double qualToErrorProbLog10(byte qual) {
        return qualToErrorProbLog10Cache[(int) qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    public static int fastRound(double d) {
        return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
    }

    public static double approximateLog10SumLog10(double a, double b, double c) {
        return approximateLog10SumLog10(a, approximateLog10SumLog10(b, c));
    }

    public static double approximateLog10SumLog10(double small, double big) {
        // make sure small is really the smaller value
        if (small > big) {
            final double t = big;
            big = small;
            small = t;
        }

        if (small == Double.NEGATIVE_INFINITY || big == Double.NEGATIVE_INFINITY)
            return big;

        final double diff = big - small;
        if (diff >= MAX_JACOBIAN_TOLERANCE)
            return big;

        // OK, so |y-x| < tol: we use the following identity then:
        // we need to compute log10(10^x + 10^y)
        // By Jacobian logarithm identity, this is equal to
        // max(x,y) + log10(1+10^-abs(x-y))
        // we compute the second term as a table lookup with integer quantization
        // we have pre-stored correction for 0,0.1,0.2,... 10.0
        final int ind = fastRound(diff * JACOBIAN_LOG_TABLE_INV_STEP); // hard rounding
        return big + jacobianLogTable[ind];
    }

    public double computeReadLikelihoodGivenHaplotype(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals,
                                                      final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = readBases.length + 2;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 2;

        // initial arrays to hold the probabilities of being in the match, insertion and deletion cases
        final double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        initializeArrays(matchMetricArray, XMetricArray, YMetricArray, X_METRIC_LENGTH);

        return computeReadLikelihoodGivenHaplotype(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, 0, matchMetricArray, XMetricArray, YMetricArray);
    }

    public double computeReadLikelihoodGivenHaplotype(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals,
                                                      final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, final int hapStartIndex,
                                                      final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = readBases.length + 2;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 2;

        // ensure that all the qual scores have valid values
        for (int iii = 0; iii < readQuals.length; iii++) {
            readQuals[iii] = (readQuals[iii] < MIN_USABLE_Q_SCORE ? MIN_USABLE_Q_SCORE : (readQuals[iii] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[iii]));
        }

		// simple rectangular version of update loop, slow
		for (int iii = 1; iii < X_METRIC_LENGTH; iii++) {
			for (int jjj = hapStartIndex + 1; jjj < Y_METRIC_LENGTH; jjj++) {
				if ((iii == 1 && jjj == 1)) {
					continue;
				}
				updateCell(iii, jjj, haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP,
						matchMetricArray, XMetricArray, YMetricArray);
			}
		}

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = X_METRIC_LENGTH - 1;
        final int endJ = Y_METRIC_LENGTH - 1;
        return approximateLog10SumLog10(matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]);
    }

    private void updateCell(final int indI, final int indJ, final byte[] haplotypeBases, final byte[] readBases,
                            final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP,
                            final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray) {

        // the read and haplotype indices are offset by one because the state arrays have an extra column to hold the initial conditions
        final int im1 = indI - 1;
        final int jm1 = indJ - 1;

        // update the match array
        double pBaseReadLog10 = 0.0; // Math.log10(1.0);
        if (im1 > 0 && jm1 > 0) { // the emission probability is applied when leaving the state
            final byte x = readBases[im1 - 1];
            final byte y = haplotypeBases[jm1 - 1];
            final byte qual = readQuals[im1 - 1];
            pBaseReadLog10 = (x == y || x == (byte) 'N' || y == (byte) 'N' ? qualToProbLog10(qual) : qualToErrorProbLog10(qual));
        }
        final int qualIndexGOP = (im1 == 0 ? DEFAULT_GOP + DEFAULT_GOP : (insertionGOP[im1 - 1] + deletionGOP[im1 - 1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : insertionGOP[im1 - 1] + deletionGOP[im1 - 1]));
        final double d0 = qualToProbLog10((byte) qualIndexGOP);
        final double e0 = (im1 == 0 ? qualToProbLog10(DEFAULT_GCP) : qualToProbLog10(overallGCP[im1 - 1]));
        matchMetricArray[indI][indJ] = pBaseReadLog10 + approximateLog10SumLog10(matchMetricArray[indI - 1][indJ - 1] + d0, XMetricArray[indI - 1][indJ - 1] + e0, YMetricArray[indI - 1][indJ - 1] + e0);

        // update the X (insertion) array
        final double d1 = (im1 == 0 ? qualToErrorProbLog10(DEFAULT_GOP) : qualToErrorProbLog10(insertionGOP[im1 - 1]));
        final double e1 = (im1 == 0 ? qualToErrorProbLog10(DEFAULT_GCP) : qualToErrorProbLog10(overallGCP[im1 - 1]));
        final double qBaseReadLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        XMetricArray[indI][indJ] = qBaseReadLog10 + approximateLog10SumLog10(matchMetricArray[indI - 1][indJ] + d1, XMetricArray[indI - 1][indJ] + e1);

        // update the Y (deletion) array, with penalty of zero on the left and right flanks to allow for a local alignment within the haplotype
        final double d2 = (im1 == 0 || im1 == readBases.length ? 0.0 : qualToErrorProbLog10(deletionGOP[im1 - 1]));
        final double e2 = (im1 == 0 || im1 == readBases.length ? 0.0 : qualToErrorProbLog10(overallGCP[im1 - 1]));
        final double qBaseRefLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        YMetricArray[indI][indJ] = qBaseRefLog10 + approximateLog10SumLog10(matchMetricArray[indI][indJ - 1] + d2, YMetricArray[indI][indJ - 1] + e2);
    }

    /*************************************************************************
     * Classes and Functions used to parse and apply the tests to the
     * PairHMM class
     ************************************************************************/

    private static class XReadLines implements Iterator<String>, Iterable<String> {
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

    private static void runTests(LinkedList<Map<String, byte[]>> testCache, PairHMM hmm, boolean debug, boolean write) {
        if (write) {
            try {

                List<Double> results = new LinkedList<Double>();

                for (Map<String, byte[]> test : testCache) {

                    double result = hmm.computeReadLikelihoodGivenHaplotype(test.get("reference"), test.get("read"),
                            test.get("baseQuals"), test.get("insQuals"), test.get("delQuals"), test.get("gcps"));

                    if (debug) {
                        System.out.printf(" Result:%4.2f%n" +
                                "==========================================================%n", result);
                    }
                    results.add(result);
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output.txt"), true));
                for (Double result : results) {
                    bw.write(String.format("%f", result));
                    bw.newLine();
                }
                bw.close();

                System.out.printf("%d - First run-through complete.%n", System.currentTimeMillis() - startTime);
            } catch (Exception e) {//Catch exception if any
                throw new RuntimeException(e);
            }
        } else {
            for (Map<String, byte[]> test : testCache) {

                double result = hmm.computeReadLikelihoodGivenHaplotype(test.get("reference"), test.get("read"),
                        test.get("baseQuals"), test.get("insQuals"), test.get("delQuals"), test.get("gcps"));

                if (debug) {
                    System.out.printf(" Result:%4.2f%n" +
                            "==========================================================%n", result);
                }
            }
            System.out.printf("%d - No I/O run-through complete.%n", System.currentTimeMillis() - startTime);
        }
    }

    private static LinkedList<Map<String, byte[]>> parseFile(String filePath, boolean debug) {
        LinkedList<Map<String, byte[]>> testCache = new LinkedList<Map<String, byte[]>>();
        try {

            for (String line : new PairHMM.XReadLines(new File(filePath))) {
                String[] p = line.split(" ");

                if (debug) {
                    System.out.println("REF:" + p[0]);
                    System.out.println("MUT:" + p[1]);
                    System.out.println("BQ: " + p[2]);
                    System.out.println("IQ: " + p[3]);
                    System.out.println("DQ: " + p[4]);
                    System.out.println("GCP:" + p[5]);
                }

                Map<String, byte[]> test = new HashMap<String, byte[]>();
                test.put("reference", p[0].getBytes());
                test.put("read", p[1].getBytes());
                test.put("baseQuals", getQuals(p[2]));
                test.put("insQuals", getQuals(p[3]));
                test.put("delQuals", getQuals(p[4]));
                test.put("gcps", getQuals(p[5]));

                testCache.add(test);
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        return testCache;
    }

    static byte[] getQuals(String quals) {
        byte[] output = quals.getBytes();

        for (int i = 0; i < output.length; i++) {
            output[i] -= (byte) ' ';
        }

        return output;
    }

    public static void main(String[] argv) {
        boolean debug = false;
        PairHMM hmm = new PairHMM();

        List<String> args = new LinkedList<String>(Arrays.asList(argv));
        if (args.contains("-debug")) {
            debug = true;
            args.remove("-debug");
        }
        if (args.size() < 1) {
            throw new RuntimeException("\r\nYou must specify a file name for input.Format below:\n" +
                    "hap ref basequals insquals delquals gcp\n ----------------------------\n" +
                    "Run with -debug for debug output");
        } else {
            System.out.printf("%d - Loading input data ...%n", System.currentTimeMillis() - startTime);
            LinkedList<Map<String, byte[]>> testCache = new LinkedList<Map<String, byte[]>>();
            for (String arg : args) {
                testCache.addAll(parseFile(arg, debug));
            }
            System.out.printf("%d - DONE loading input data! Now calculating%n", System.currentTimeMillis() - startTime);
            boolean firstTime = true;
            for (int i = 0; i < 100; i++) {
                runTests(testCache, hmm, debug, firstTime);
                firstTime = false;
            }

        }

    }




}
