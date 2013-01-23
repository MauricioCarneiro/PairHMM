import java.io.BufferedWriter;
import java.io.File;
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
    private static LinkedList<TestRow> testsPerHaplotype = new LinkedList<TestRow>();
    private static double PRECISION = 0.001;

    public RunTest() {
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

        return pairHMM.computeReadLikelihoodGivenHaplotypeLog10( haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);
    }


    private static boolean runTestsHelper(List<Double> results, boolean debug) {

        for (TestRow test : testsPerHaplotype) {
            Double result = hmm(test.getHaplotypeBases(), test.getReadBases(),
                test.getReadQuals(), test.getReadInsQuals(),
                test.getReadDelQuals(), test.getOverallGCP(),
                test.getHaplotypeStart(), test.getReachedReadValye());

            if (debug) {
                System.out.printf(" Result:%4.3f%n" +
                        "==========================================================%n", result);
            }
            if (test.getLikelihood() != result) {
                System.out.println("Wrong result. Expected " + test.getLikelihood() + " , actual: " + result);
                return false;
            }
            results.add(result);
        }

        // clear = new haplotype will start
        testsPerHaplotype.clear();
        return true;
    }

    private static boolean runTests(LinkedList<Map<String, TestRow>> testCache, boolean debug, boolean write) {

        if (write) {
            try {

                List<Double> results = new LinkedList<Double>();
                byte[] currentHaplotype = {};
                int X_METRIC_LENGTH = 0;
                int Y_METRIC_LENGTH;

                for (Map<String, TestRow> test : testCache) {
                    TestRow currentTest = test.get("testInstance");

                    if (X_METRIC_LENGTH < currentTest.getReadBases().length)   {
                        X_METRIC_LENGTH = currentTest.getReadBases().length;
                    }
                    // first pass
                    if(currentHaplotype.length == 0) {
                        currentHaplotype = currentTest.getHaplotypeBases();
                    }

                    // check if new haplotype
                    if(currentTest.getHaplotypeStart() == 0 & !Arrays.equals(currentHaplotype, currentTest.getHaplotypeBases()) ) {
                        Y_METRIC_LENGTH = currentHaplotype.length;
                        pairHMM.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);

                        runTestsHelper(results, debug);

                        currentHaplotype = currentTest.getHaplotypeBases();
                    }
                    testsPerHaplotype.add(currentTest);
                }

                Y_METRIC_LENGTH = currentHaplotype.length;
                pairHMM.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);
                runTestsHelper(results, debug);

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output.txt"), false));
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
            for (Map<String, TestRow> test : testCache) {
                TestRow currentTest = test.get("testInstance");
                double result = hmm(currentTest.getHaplotypeBases(), currentTest.getReadBases(),
                        currentTest.getReadQuals(), currentTest.getReadInsQuals(),
                        currentTest.getReadDelQuals(), currentTest.getOverallGCP(),
                        currentTest.getHaplotypeStart(), currentTest.getReachedReadValye());

                if ((currentTest.getLikelihood() - result) < PRECISION) {
                    System.out.println("Wrong result. Expected " + currentTest.getLikelihood() + " , actual: " + result);
                    return false;
                }

                if (debug) {
                    System.out.printf(" Result:%4.2f%n" +
                            "==========================================================%n", result);
                }
            }
            System.out.printf("%d - No I/O run-through complete.%n", System.currentTimeMillis() - startTime);
        }
        return true;
    }

    private static LinkedList<Map<String, TestRow>> parseFile(String filePath, boolean debug) {
        LinkedList<Map<String, TestRow>> testCache = new LinkedList<Map<String, TestRow>>();
        try {
            for (String line : new XReadLines(new File(filePath))) {
                String[] p = line.split(" ");

                if (debug) {
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
                boolean reachedRead;
                if (p[7].trim().equals("true")) {
                    reachedRead = true;
                } else {
                    reachedRead = false;
                }

                Map<String, TestRow> test = new HashMap<String, TestRow>();
                TestRow testRow = new TestRow(p[0].getBytes(), p[1].getBytes(),
                        convertToByteArray(p[2]), convertToByteArray(p[3]),
                        convertToByteArray(p[4]), convertToByteArray(p[5]),
                        Integer.parseInt(p[6]), reachedRead,
                        Double.parseDouble(p[8]));

                test.put("testInstance", testRow);
                testCache.add(test);
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        return testCache;
    }

    static byte[] convertToByteArray(String quals) {
        byte[] output = quals.getBytes();

        for (int i = 0; i < output.length; i++) {
            output[i] -= 33;
        }

        return output;
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
                    "hap ref basequals insquals delquals gcp\n ----------------------------\n" +
                    "Run with -debug for debug output");
        } else {
            System.out.printf("%d - Loading input data ...%n", System.currentTimeMillis() - startTime);
            LinkedList<Map<String, TestRow>> testCache = new LinkedList<Map<String, TestRow>>();
            for (String arg : args) {
                testCache.addAll(parseFile(arg, debug));
            }
            System.out.printf("%d - DONE loading input data! Now calculating%n", System.currentTimeMillis() - startTime);
            if(runTests(testCache, debug, true)) {
                System.out.println("Tests successful");
            } else {
                System.out.println("Tests unsuccessful");
            }
        }

    }
}
