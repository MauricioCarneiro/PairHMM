
/**
 * @author Mauricio Carneiro
 * @since 8/20/12
 */
public class QualityUtils {
    private static final int MAXN = 50000;
    private static double[] qualToErrorProbLog10Cache = new double[256];
    private static double[] qualToProbLog10Cache = new double[256];
    private static final double[] log10Cache;
    private static final double[] log10FactorialCache;
    private static final int LOG10_CACHE_SIZE = 4 * MAXN;  // we need to be able to go up to 2*(2N) when calculating some of the coefficients
    public final static byte MIN_USABLE_Q_SCORE = 6;

    static {
        for (int i = 0; i < 256; i++)
            qualToErrorProbLog10Cache[i] = qualToErrorProbLog10Raw(i);
        for (int i = 0; i < 256; i++)
            qualToProbLog10Cache[i] = qualToProbLog10Raw(i);

        log10Cache = new double[LOG10_CACHE_SIZE];
        log10FactorialCache = new double[LOG10_CACHE_SIZE];

        log10Cache[0] = Double.NEGATIVE_INFINITY;
        for (int k = 1; k < LOG10_CACHE_SIZE; k++) {
            log10Cache[k] = Math.log10(k);
            log10FactorialCache[k] = log10FactorialCache[k - 1] + log10Cache[k];
        }

    }

    /**
     * Convert a quality score to a probability.  This is the Phred-style
     * conversion, *not* the Illumina-style conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0-255)
     * @return a probability (0.0-1.0)
     */
    static public double qualToProb(byte qual) {
        return 1.0 - qualToErrorProb(qual);
    }

    static private double qualToProbLog10Raw(int qual) {
        return Math.log10(1.0 - qualToErrorProbRaw(qual));
    }

    static public double qualToProbLog10(byte qual) {
        return qualToProbLog10Cache[(int) qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    /**
     * Convert a quality score to a probability of error.  This is the Phred-style conversion, *not* the Illumina-style
     * conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0 - 255)
     * @return a probability (0.0 - 1.0)
     */
    static private double qualToErrorProbRaw(int qual) {
        return qualToErrorProb((double) qual);
    }

    public static double qualToErrorProb(final double qual) {
        return Math.pow(10.0, (qual) / -10.0);
    }

    static private double qualToErrorProbLog10Raw(int qual) {
        return ((double) qual) / -10.0;
    }

    static public double qualToErrorProbLog10(byte qual) {
        return qualToErrorProbLog10Cache[(int) qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }
}
