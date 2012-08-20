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
    private final static byte MIN_USABLE_Q_SCORE = 6;
    private static final double [] firstRowConstantMatrix = {
            QualityUtils.qualToProbLog10((byte) (DEFAULT_GOP + DEFAULT_GOP)),
            QualityUtils.qualToProbLog10(DEFAULT_GCP),
            QualityUtils.qualToErrorProbLog10(DEFAULT_GOP),
            QualityUtils.qualToErrorProbLog10(DEFAULT_GCP),
            0.0,
            0.0
    };

    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     * @param distanceMatrix the pre-allocated distance matrix
     */
    public static void initializeDistanceMatrix(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, int startIndex, double [][] distanceMatrix) {
        // initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                distanceMatrix[i+2][j+2] = (x == y || x == (byte) 'N' || y == (byte) 'N' ? QualityUtils.qualToProbLog10(qual) : QualityUtils.qualToErrorProbLog10(qual));
            }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     * @param constantMatrix the pre-allocated constant matrix
     */
    public static void initializeConstants(final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, double [][] constantMatrix) {
        final int l = insertionGOP.length;
        constantMatrix[1] = firstRowConstantMatrix;
        for (int i = 0; i < l; i++) {
            final int qualIndexGOP = Math.min(insertionGOP[i] + deletionGOP[i], MAX_CACHED_QUAL);
            constantMatrix[i+2][0] = QualityUtils.qualToProbLog10((byte) qualIndexGOP);
            constantMatrix[i+2][1] = QualityUtils.qualToProbLog10(overallGCP[i]);
            constantMatrix[i+2][2] = QualityUtils.qualToErrorProbLog10(insertionGOP[i]);
            constantMatrix[i+2][3] = QualityUtils.qualToErrorProbLog10(overallGCP[i]);
            constantMatrix[i+2][4] = QualityUtils.qualToErrorProbLog10(deletionGOP[i]);
            constantMatrix[i+2][5] = QualityUtils.qualToErrorProbLog10(overallGCP[i]);
        }
        constantMatrix[l+1][4] = 0.0;
        constantMatrix[l+1][5] = 0.0;
    }

    public static void initializeArrays(final int readDimension, final int hapDimension, final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray) {

        // fill the matrix with -inf
        for (int i = 0; i < readDimension; i++) {
            Arrays.fill(matchMetricArray[i], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[i], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[i], Double.NEGATIVE_INFINITY);
        }
        // the initial condition
        matchMetricArray[1][1] = 0.0; // Math.log10(1.0);

        // fill in the first row
        for (int j = 2; j < hapDimension; j++) {
            updateCell(1, j, 0.0, firstRowConstantMatrix, matchMetricArray, XMetricArray, YMetricArray);
        }
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
    public static double computeReadLikelihoodGivenHaplotype(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        // ensure that all the qual scores have valid values
        for (int i = 0; i < readQuals.length; i++) {
            readQuals[i] = (readQuals[i] < MIN_USABLE_Q_SCORE ? MIN_USABLE_Q_SCORE : (readQuals[i] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[i]));
        }

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int readDimension = readBases.length + 2;
        final int haplotypeDimension = haplotypeBases.length + 2;

        // initial arrays to hold the probabilities of being in the match, insertion and deletion cases
        final double[][] matchMetricArray = new double[readDimension][haplotypeDimension];
        final double[][] XMetricArray = new double[readDimension][haplotypeDimension];
        final double[][] YMetricArray = new double[readDimension][haplotypeDimension];
        final double[][] constantMatrix = new double[readDimension][6];
        final double[][] distanceMatrix = new double[readDimension][haplotypeDimension];
        initializeDistanceMatrix(haplotypeBases, readBases, readQuals, 0, distanceMatrix);
        initializeConstants(insertionGOP, deletionGOP, overallGCP, constantMatrix);
        initializeArrays(readDimension, haplotypeDimension, matchMetricArray, XMetricArray, YMetricArray);
        return computeReadLikelihoodGivenHaplotype( readDimension, haplotypeDimension, 0, matchMetricArray, XMetricArray, YMetricArray, constantMatrix, distanceMatrix);
    }

    /**
     * Computes the Pair HMM matrix for a previously initialized matrix.
     *
     * Only use this method with startIndex >= 2, if you need to fill in the first two rows of the matrix, use the
     * other version that initializes and computes the matrix.
     *
     * This startIndex is equivalent to the haplotype start index plus one because of the padding in the left of the matrices
     *   startIndex = hapStartIndex + 1
     *
     * (startIndex can never be < 2, indices below 2 are initialized in the initializeArrays routine)
     *
     * @param readMetricLength length of the read dimenstion of the metrics arrays (readBases + 2)
     * @param hapMetricLength  length of the haplotype dimension of the metrics arrays (haplotypeBases + 2)
     * @param startIndex       haplotype index to start the calculation (assumes all previous indices have already been filled on the three metric arrays)
     * @param matchMetricArray partially pre-filled matches metric array
     * @param XMetricArray     partially pre-filled insertion metric array
     * @param YMetricArray     partially pre-filled deletion metric array
     * @param constantMatrix   the matrix containing all the read constants
     * @param distanceMatrix   probability constants for all combinations of i and j between haplotype and read
     * @return the likelihood of the alignment between read and haplotype
     */
    public static double computeReadLikelihoodGivenHaplotype(int readMetricLength, int hapMetricLength, final int startIndex, final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray, double[][] constantMatrix, double[][] distanceMatrix) {
        for (int i = 2; i < readMetricLength; i++) {
            for (int j = startIndex+1; j < hapMetricLength; j++) {
                updateCell(i, j, distanceMatrix[i][j], constantMatrix[i], matchMetricArray, XMetricArray, YMetricArray);
            }
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = readMetricLength - 1;
        final int endJ = hapMetricLength - 1;
        return MathUtils.approximateLog10SumLog10(matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]);
    }

    /**
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions

     * @param indI             row index in the matrices to update
     * @param indJ             column index in the matrices to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param constants        an array with the six constants relevant to this location
     * @param matchMetricArray the matches likelihood matrix
     * @param XMetricArray     the insertions likelihood matrix
     * @param YMetricArray     the deletions likelihood matrix
     */
    private static void updateCell(final int indI, final int indJ, double prior, double[] constants, final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray) {
        matchMetricArray[indI][indJ] = prior + MathUtils.approximateLog10SumLog10(matchMetricArray[indI - 1][indJ - 1] + constants[0], XMetricArray[indI - 1][indJ - 1] + constants[1], YMetricArray[indI - 1][indJ - 1] + constants[1]);
        XMetricArray[indI][indJ] = MathUtils.approximateLog10SumLog10(matchMetricArray[indI - 1][indJ] + constants[2], XMetricArray[indI - 1][indJ] + constants[3]);
        YMetricArray[indI][indJ] = MathUtils.approximateLog10SumLog10(matchMetricArray[indI][indJ - 1] + constants[4], YMetricArray[indI][indJ - 1] + constants[5]);
    }

}
