package org.broadinstitute.utils;

/**
 * @author Mauricio Carneiro
 * @since 8/16/12
 */
public class MathUtils {

    private static final float[] jacobianLogTable;
    private static final float JACOBIAN_LOG_TABLE_STEP = 0.001f;
    private static final float MAX_JACOBIAN_TOLERANCE = 8.f;
    private static final float JACOBIAN_LOG_TABLE_INV_STEP = 1.f / 0.001f;
    private static final int JACOBIAN_LOG_TABLE_SIZE = (int) (MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1;

    static {
        jacobianLogTable = new float[JACOBIAN_LOG_TABLE_SIZE];

        for (int k = 0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
            jacobianLogTable[k] = (float) Math.log10(1.0 + Math.pow(10.0, -((float) k) * JACOBIAN_LOG_TABLE_STEP));

        }
    }

    /**
     * Checks that the result is a well-formed log10 probability
     *
     * @param result a supposedly well-formed log10 probability value.  By default allows
     *               -Infinity values, as log10(0.0) == -Infinity.
     * @return true if result is really well formed
     */
    public static boolean goodLog10Probability(final float result) {
        return goodLog10Probability(result, true);
    }

    /**
     * Checks that the result is a well-formed log10 probability
     *
     * @param result a supposedly well-formed log10 probability value
     * @param allowNegativeInfinity should we consider a -Infinity value ok?
     * @return true if result is really well formed
     */
    public static boolean goodLog10Probability(final float result, final boolean allowNegativeInfinity) {
        return result <= 0.0 && result != Float.POSITIVE_INFINITY && (allowNegativeInfinity || result != Float.NEGATIVE_INFINITY) && ! Float.isNaN(result);
    }

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    public static int fastRound(float d) {
        return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
    }

    public static int maxElementIndex(final float[] array, final int endIndex) {
        if (array == null || array.length == 0)
            throw new IllegalArgumentException("Array cannot be null!");

        int maxI = 0;
        for (int i = 1; i < endIndex; i++) {
            if (array[i] > array[maxI])
                maxI = i;
        }

        return maxI;
    }

    public static int maxElementIndex(final float[] array) {
        return maxElementIndex(array, array.length);
    }

    public static float arrayMax(final float[] array, final int endIndex) {
        return array[maxElementIndex(array, endIndex)];
    }

    public static float log10sumLog10(float[] log10p, int start, int finish) {
        float sum = 0.f;

        float maxValue = arrayMax(log10p, finish);
        if(maxValue == Float.NEGATIVE_INFINITY)
            return maxValue;

        for (int i = start; i < finish; i++) {
            sum += (float) Math.pow(10.0, log10p[i] - maxValue);
        }

        return (float) Math.log10(sum) + maxValue;
    }

    public static float approximateLog10SumLog10(final float[] vals) {
        return approximateLog10SumLog10(vals, vals.length);
    }

    public static float approximateLog10SumLog10(final float[] vals, final int endIndex) {

        final int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
        float approxSum = vals[maxElementIndex];

        for (int i = 0; i < endIndex; i++) {
            if (i == maxElementIndex || vals[i] == Float.NEGATIVE_INFINITY)
                continue;

            final float diff = approxSum - vals[i];
            if (diff < MathUtils.MAX_JACOBIAN_TOLERANCE) {
                // See notes from the 2-inout implementation below
                final int ind = fastRound(diff * MathUtils.JACOBIAN_LOG_TABLE_INV_STEP); // hard rounding
                approxSum += MathUtils.jacobianLogTable[ind];
            }
        }

        return approxSum;
    }

    public static float log10sumLog10(float[] log10p, int start) {
        return log10sumLog10(log10p, start, log10p.length);
    }

    public static float log10sumLog10(float[] log10values) {
        return log10sumLog10(log10values, 0);
    }

    public static float approximateLog10SumLog10(float small, float big) {
        // make sure small is really the smaller value
        if (small > big) {
            final float t = big;
            big = small;
            small = t;
        }

        if (small == Float.NEGATIVE_INFINITY)
            return big;

        final float diff = big - small;
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
}
