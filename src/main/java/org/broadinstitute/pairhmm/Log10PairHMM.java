/*
 * Copyright (c) 2012 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.pairhmm;

import org.broadinstitute.utils.MathUtils;
import org.broadinstitute.utils.QualityUtils;

import java.util.Arrays;

/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 *
 * User: rpoplin
 * Date: 3/1/12
 */
public final class Log10PairHMM extends PairHMM {
    /**
     * Should we use exact log10 calculation (true), or an approximation (false)?
     */
    private final boolean doExactLog10;
    /**
     * Create an uninitialized PairHMM
     *
     * @param doExactLog10 should the log10 calculations be exact (slow) or approximate (faster)
     */
    public Log10PairHMM(final boolean doExactLog10) {
        this.doExactLog10 = doExactLog10;
    }

    /**
     * Is this HMM using exact log10 calculations?
     * @return true if exact, false if approximate
     */
    public boolean isDoingExactLog10Calculations() {
        return doExactLog10;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize(final int readMaxLength,  final int haplotypeMaxLength) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        for( int iii=0; iii < paddedMaxReadLength; iii++ ) {
            Arrays.fill(matchMatrix[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(insertionMatrix[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(deletionMatrix[iii], Double.NEGATIVE_INFINITY);
        }

        transition = new double[paddedMaxReadLength][6];
        prior = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               final int hapStartIndex,
                                                               final boolean recacheReadValues ) {

        // don't have to initialize every time... revise this!
        initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);
        initializeProbabilities(Math.log10(1.0/haplotypeBases.length), insertionGOP, deletionGOP, overallGCP);

        // simple rectangular version of update loop, slow
        for( int i = 1; i < paddedReadLength; i++ ) {
            for( int j = hapStartIndex + 1; j < paddedHaplotypeLength; j++ ) {
                updateCell(i, j, prior[i][j], transition[i]);
            }
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = paddedReadLength - 1;
        double result = myLog10SumLog10(new double[]{matchMatrix[endI][1], insertionMatrix[endI][1]});
        for (int j = 2; j < paddedHaplotypeLength; j++)
            result = myLog10SumLog10(new double[]{result, matchMatrix[endI][j], insertionMatrix[endI][j]});

        return result;
    }

    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     */
    public void initializePriors(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {
        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProbLog10(qual) : QualityUtils.qualToErrorProbLog10(qual) );
            }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    private void initializeProbabilities(final double initialValue, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
        for (int j = 0; j < paddedHaplotypeLength; j++) {
            deletionMatrix[0][j] = initialValue;
        }

        for (int i = 0; i < insertionGOP.length; i++) {
            final int qualIndexGOP = Math.min(insertionGOP[i] + deletionGOP[i], Byte.MAX_VALUE);
            transition[i+1][0] = QualityUtils.qualToProbLog10((byte) qualIndexGOP);
            transition[i+1][1] = QualityUtils.qualToProbLog10(overallGCP[i]);
            transition[i+1][2] = QualityUtils.qualToErrorProbLog10(insertionGOP[i]);
            transition[i+1][3] = QualityUtils.qualToErrorProbLog10(overallGCP[i]);
            transition[i+1][4] = QualityUtils.qualToErrorProbLog10(deletionGOP[i]);
            transition[i+1][5] = QualityUtils.qualToErrorProbLog10(overallGCP[i]);
        }

        // note that we initialized the constants
        constantsAreInitialized = true;
    }


    /**
     * Compute the log10SumLog10 of the values
     *
     * NOTE NOTE NOTE
     *
     * Log10PairHMM depends critically on this function tolerating values that are all -Infinity
     * and the sum returning -Infinity.  Note good.  Needs to be fixed.
     *
     * NOTE NOTE NOTE
     *
     * @param values an array of log10 probabilities that need to be summed
     * @return the log10 of the sum of the probabilities
     */
    private double myLog10SumLog10(final double[] values) {
        return doExactLog10 ? MathUtils.log10sumLog10(values) : MathUtils.approximateLog10SumLog10(values);
    }

    /**
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions

     * @param indI             row index in the matrices to update
     * @param indJ             column index in the matrices to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param transitition        an array with the six transitition relevant to this location
     */
    private void updateCell( final int indI, final int indJ, final double prior, final double[] transitition) {

        matchMatrix[indI][indJ] = prior +
                myLog10SumLog10(new double[]{matchMatrix[indI - 1][indJ - 1] + transitition[0],
                        insertionMatrix[indI - 1][indJ - 1] + transitition[1],
                        deletionMatrix[indI - 1][indJ - 1] + transitition[1]});
        insertionMatrix[indI][indJ] = myLog10SumLog10(new double[] {matchMatrix[indI - 1][indJ] + transitition[2], insertionMatrix[indI - 1][indJ] + transitition[3]});
        deletionMatrix[indI][indJ]  = myLog10SumLog10(new double[] {matchMatrix[indI][indJ - 1] + transitition[4],  deletionMatrix[indI][indJ - 1] + transitition[5]});
    }
}
