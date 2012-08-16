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

            for (String line : new XReadLines(new File(filePath))) {
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
            runTests(testCache, hmm, debug, true);
        }

    }
}
