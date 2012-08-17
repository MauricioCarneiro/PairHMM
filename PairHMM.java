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

import java.util.Arrays;

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

    public static void initializeArrays(final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray, final int X_METRIC_LENGTH, final int Y_METRIC_LENGTH) {

        for (int i = 0; i < X_METRIC_LENGTH; i++) {
            Arrays.fill(matchMetricArray[i], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[i], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[i], Double.NEGATIVE_INFINITY);
        }
        // the initial condition
        matchMetricArray[1][1] = 0.0; // Math.log10(1.0);

        final double[] d = {MathUtils.qualToProbLog10((byte) (DEFAULT_GOP + DEFAULT_GOP)), MathUtils.qualToErrorProbLog10(DEFAULT_GOP), 0.0};
        final double[] e = {MathUtils.qualToProbLog10(DEFAULT_GCP), MathUtils.qualToErrorProbLog10(DEFAULT_GCP), 0.0};

        for (int j = 2; j < Y_METRIC_LENGTH; j++) {
            updateCell(1, j, 0.0, d, e, matchMetricArray, XMetricArray, YMetricArray);
        }
    }

    public static double finalizeArrays(final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, int hapStartIndex, final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int X_METRIC_LENGTH, final int Y_METRIC_LENGTH) {
        final int i = X_METRIC_LENGTH - 1;
        final int qualIndexGOP = Math.min(insertionGOP[i - 2] + deletionGOP[i - 2], MAX_CACHED_QUAL);
        final double d[] = new double[3];
        final double e[] = new double[3];
        d[0] = MathUtils.qualToProbLog10((byte) qualIndexGOP);
        e[0] = MathUtils.qualToProbLog10(overallGCP[i - 2]);
        d[1] = MathUtils.qualToErrorProbLog10(insertionGOP[i - 2]);
        e[1] = MathUtils.qualToErrorProbLog10(overallGCP[i - 2]);
        d[2] = 0.0;
        e[2] = 0.0;
        final byte x = readBases[i - 2];
        final byte qual = readQuals[i - 2];
        for (int j = Math.max(2, hapStartIndex + 1); j < Y_METRIC_LENGTH; j++) {
            final byte y = haplotypeBases[j - 2];
            final double pBaseReadLog10 = (x == y || x == (byte) 'N' || y == (byte) 'N' ? MathUtils.qualToProbLog10(qual) : MathUtils.qualToErrorProbLog10(qual));
            updateCell(i, j, pBaseReadLog10, d, e, matchMetricArray, XMetricArray, YMetricArray);
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = X_METRIC_LENGTH - 1;
        final int endJ = Y_METRIC_LENGTH - 1;
        return MathUtils.approximateLog10SumLog10(matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]);

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
    public double computeReadLikelihoodGivenHaplotype(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {

        // ensure that all the qual scores have valid values
        for (int i = 0; i < readQuals.length; i++) {
            readQuals[i] = (readQuals[i] < MIN_USABLE_Q_SCORE ? MIN_USABLE_Q_SCORE : (readQuals[i] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[i]));
        }

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = readBases.length + 2;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 2;

        // initial arrays to hold the probabilities of being in the match, insertion and deletion cases
        final double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        initializeArrays(matchMetricArray, XMetricArray, YMetricArray, X_METRIC_LENGTH, Y_METRIC_LENGTH);
        return computeReadLikelihoodGivenHaplotype(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, 2, matchMetricArray, XMetricArray, YMetricArray);
    }

    /**
     * Computes the Pair HMM matrix for a previously initialized matrix.
     *
     * Only use this method with hapStartIndex >= 2, if you need to fill in the first two rows of the matrix, use the
     * other version that initializes and computes the matrix.
     *
     * @param haplotypeBases   reference sequence bases
     * @param readBases        comparison haplotype bases
     * @param readQuals        comparison haplotype base quals (phred-scaled)
     * @param insertionGOP     comparison haplotype insertion quals (phred-scaled)
     * @param deletionGOP      comparison haplotype deletion quals (phred-scaled)
     * @param overallGCP       comparison haplotype gap continuation quals (phred-scaled)
     * @param hapStartIndex    haplotype index to start the calculation (assumes all previous indices have already been
     *                         filled on the three metric arrays)
     * @param matchMetricArray partially pre-filled matches metric array
     * @param XMetricArray     partially pre-filled insertion metric array
     * @param YMetricArray     partially pre-filled deletion metric array
     * @return the likelihood of the alignment between read and haplotype
     */
    public double computeReadLikelihoodGivenHaplotype(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, final int hapStartIndex, final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = readBases.length + 2;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 2;

        // simple rectangular version of update loop, slow
        for (int i = hapStartIndex; i < X_METRIC_LENGTH - 1; i++) {

            final int qualIndexGOP = Math.min(insertionGOP[i - 2] + deletionGOP[i - 2], MAX_CACHED_QUAL);
            final double d[] = new double[3];
            final double e[] = new double[3];
            d[0] = MathUtils.qualToProbLog10((byte) qualIndexGOP);
            e[0] = MathUtils.qualToProbLog10(overallGCP[i - 2]);
            d[1] = MathUtils.qualToErrorProbLog10(insertionGOP[i - 2]);
            e[1] = MathUtils.qualToErrorProbLog10(overallGCP[i - 2]);
            d[2] = MathUtils.qualToErrorProbLog10(deletionGOP[i - 2]);
            e[2] = MathUtils.qualToErrorProbLog10(overallGCP[i - 2]);

            // the emission probability is applied when leaving the state
            final byte x = readBases[i - 2];
            final byte qual = readQuals[i - 2];

            // In case hapStart > 0, we will unnecessarily call this method (avoiding an if statement)
            updateCell(i, 1, 0.0, d, e, matchMetricArray, XMetricArray, YMetricArray);

            for (int j = hapStartIndex; j < Y_METRIC_LENGTH; j++) {
                final byte y = haplotypeBases[j - 2];
                final double pBaseReadLog10 = (x == y || x == (byte) 'N' || y == (byte) 'N' ? MathUtils.qualToProbLog10(qual) : MathUtils.qualToErrorProbLog10(qual));
                updateCell(i, j, pBaseReadLog10, d, e, matchMetricArray, XMetricArray, YMetricArray);
            }
        }

        return finalizeArrays(matchMetricArray, XMetricArray, YMetricArray, insertionGOP, deletionGOP, overallGCP, 0, haplotypeBases, readBases, readQuals, X_METRIC_LENGTH, Y_METRIC_LENGTH);
    }

    /**
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions
     */

    private static void updateCell(final int indI, final int indJ, double pBaseReadLog10, double[] d, double[] e, final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray) {
        matchMetricArray[indI][indJ] = pBaseReadLog10 + MathUtils.approximateLog10SumLog10(matchMetricArray[indI - 1][indJ - 1] + d[0], XMetricArray[indI - 1][indJ - 1] + e[0], YMetricArray[indI - 1][indJ - 1] + e[0]);
        XMetricArray[indI][indJ] = MathUtils.approximateLog10SumLog10(matchMetricArray[indI - 1][indJ] + d[1], XMetricArray[indI - 1][indJ] + e[1]);
        YMetricArray[indI][indJ] = MathUtils.approximateLog10SumLog10(matchMetricArray[indI][indJ - 1] + d[2], YMetricArray[indI][indJ - 1] + e[2]);
    }

}
