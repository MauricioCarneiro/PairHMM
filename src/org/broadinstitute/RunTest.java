
    import java.io.BufferedWriter;
    import java.io.File;
    import java.io.FileNotFoundException;
    import java.io.FileWriter;
    import java.util.*;

    /**
     * Testing routine for the PairHMM
     *
     * @author Mauricio Carneiro
     * @since 8/15/12
     */
    public class RunTest {

        private static final long startTime = System.currentTimeMillis();
        private static final PairHMM pairHMM = new LoglessCachingPairHMM();
        private static double PRECISION = 0.001;

        private static int X_METRIC_LENGTH = 10000;
        private static int Y_METRIC_LENGTH = 10000;

        public RunTest() {
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
                    if (Math.abs(currentTest.getLikelihood() - result) > PRECISION) {
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
    }
