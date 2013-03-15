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

/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 *
 * User: rpoplin
 * Date: 10/16/12
 */
public abstract class PairHMM {

    protected static final Byte MAX_CACHED_QUAL = Byte.MAX_VALUE;
    protected static final byte DEFAULT_GOP = (byte) 45;
    protected static final byte DEFAULT_GCP = (byte) 10;

    /**
     * unassign all the initializations to the arrays to allow java to garbage collect
     *
     * only clear if you are done with the PairHMM. After clearing, you will have to initialize it to reuse.
     */
    public void clear() {
        XMetricArray = null;
        YMetricArray = null;
        matchMetricArray = null;
    }

    protected double[][] matchMetricArray = null;
    protected double[][] XMetricArray = null;
    protected double[][] YMetricArray = null;
    protected int maxHaplotypeLength, maxReadLength;
    protected int X_METRIC_MAX_LENGTH, Y_METRIC_MAX_LENGTH;
    protected boolean initialized = false;

    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     */
    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) {
        if ( readMaxLength <= 0 ) throw new IllegalArgumentException("READ_MAX_LENGTH must be > 0 but got " + readMaxLength);
        if ( haplotypeMaxLength <= 0 ) throw new IllegalArgumentException("HAPLOTYPE_MAX_LENGTH must be > 0 but got " + haplotypeMaxLength);

        maxHaplotypeLength = haplotypeMaxLength;
        maxReadLength = readMaxLength;

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        X_METRIC_MAX_LENGTH = readMaxLength + 2;
        Y_METRIC_MAX_LENGTH = haplotypeMaxLength + 2;

        matchMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
        XMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
        YMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
        initialized = true;
    }

    /**
     * To be overloaded by subclasses to actually do calculation for #computeReadLikelihoodGivenHaplotypeLog10
     */
    public abstract double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                           final byte[] readBases,
                                                                           final byte[] readQuals,
                                                                           final byte[] insertionGOP,
                                                                           final byte[] deletionGOP,
                                                                           final byte[] overallGCP,
                                                                           final int hapStartIndex,
                                                                           final boolean recacheReadValues );


    /**
     * Print out the core hmm matrices for debugging
     */
    protected void dumpMatrices() {
        dumpMatrix("matchMetricArray", matchMetricArray);
        dumpMatrix("XMetricArray", XMetricArray);
        dumpMatrix("YMetricArray", YMetricArray);
    }

    /**
     * Print out in a human readable form the matrix for debugging
     * @param name the name of this matrix
     * @param matrix the matrix of values
     */
    private void dumpMatrix(final String name, final double[][] matrix) {
        System.out.printf("%s%n", name);
        for ( int i = 0; i < matrix.length; i++) {
            System.out.printf("\t%s[%d]", name, i);
            for ( int j = 0; j < matrix[i].length; j++ ) {
                if ( Double.isInfinite(matrix[i][j]) )
                    System.out.printf(" %15s", String.format("%f", matrix[i][j]));
                else
                    System.out.printf(" % 15.5e", matrix[i][j]);
            }
            System.out.println();
        }
    }
}
