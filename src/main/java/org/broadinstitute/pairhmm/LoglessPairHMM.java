package org.broadinstitute.pairhmm;

import org.broadinstitute.utils.QualityUtils;

/**
 * Features of this implementation over standard:
 *
 *  -- caching constants
 *  -- no logs throughout the routine (only the result)
 *  -- does not recalculate matrix for similar incoming haplotypes (uses hapStartIndex)
 *
 * note: this is the routine currently used by the HaplotypeCaller in the GATK.
 *
 * Created with IntelliJ IDEA.
 * User: rpoplin, carneiro
 * Date: 10/16/12
 */
public class LoglessPairHMM extends PairHMM {
    protected static final double SCALE_FACTOR_LOG10 = 300.0;

    double[][] constantMatrix = null; // The cache
    double[][] distanceMatrix = null; // The cache
    boolean constantsAreInitialized = false;

    /**
     * Cached data structure that describes the first row's edge condition in the HMM
     */
    protected static final double [] firstRowConstantMatrix = {
            QualityUtils.qualToProb((byte) (DEFAULT_GOP + DEFAULT_GOP)),
            QualityUtils.qualToProb(DEFAULT_GCP),
            QualityUtils.qualToErrorProb(DEFAULT_GOP),
            QualityUtils.qualToErrorProb(DEFAULT_GCP),
            1.0,
            1.0
    };

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        constantMatrix = new double[X_METRIC_MAX_LENGTH][6];
        distanceMatrix = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
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
        if ( ! constantsAreInitialized || recacheReadValues )
            initializeConstants( haplotypeBases.length, readBases.length, insertionGOP, deletionGOP, overallGCP );
        initializeDistanceMatrix( haplotypeBases, readBases, readQuals, hapStartIndex );

        // NOTE NOTE NOTE -- because of caching we need to only operate over X and Y according to this
        // read and haplotype lengths, not the max lengths
        final int readXMetricLength = readBases.length + 2;
        final int hapYMetricLength = haplotypeBases.length + 2;

        for (int i = 2; i < readXMetricLength; i++) {
//            dumpMatrices();
            for (int j = hapStartIndex+1; j < hapYMetricLength; j++) {
                updateCell(i, j, distanceMatrix[i][j], constantMatrix[i], matchMetricArray, XMetricArray, YMetricArray);
            }
        }
//        dumpMatrices();


        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = readXMetricLength - 1;
        final int endJ = hapYMetricLength - 1;
        return Math.log10( matchMetricArray[endI][endJ] + XMetricArray[endI][endJ] + YMetricArray[endI][endJ] ) - SCALE_FACTOR_LOG10;
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
    public void initializeDistanceMatrix( final byte[] haplotypeBases,
                                          final byte[] readBases,
                                          final byte[] readQuals,
                                          final int startIndex ) {

        // initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                distanceMatrix[i+2][j+2] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProb(qual) : QualityUtils.qualToErrorProb(qual) );
            }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param haplotypeSize the number of bases in the haplotype we are testing
     * @param readSize the number of bases in the read we are testing
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    private void initializeConstants( final int haplotypeSize,
                                      final int readSize,
                                      final byte[] insertionGOP,
                                      final byte[] deletionGOP,
                                      final byte[] overallGCP ) {
        // the initial condition -- must be here because it needs that actual read and haplotypes, not the maximum in init
        matchMetricArray[1][1] = Math.pow(10.0, SCALE_FACTOR_LOG10);

        // fill in the first row
        for( int jjj = 2; jjj < Y_METRIC_MAX_LENGTH; jjj++ ) {
            updateCell(1, jjj, 1.0, firstRowConstantMatrix, matchMetricArray, XMetricArray, YMetricArray);
        }

        final int l = insertionGOP.length;
        constantMatrix[1] = firstRowConstantMatrix;
        for (int i = 0; i < l; i++) {
            final int qualIndexGOP = Math.min(insertionGOP[i] + deletionGOP[i], Byte.MAX_VALUE);
            constantMatrix[i+2][0] = QualityUtils.qualToProb((byte) qualIndexGOP);
            constantMatrix[i+2][1] = QualityUtils.qualToProb(overallGCP[i]);
            constantMatrix[i+2][2] = QualityUtils.qualToErrorProb(insertionGOP[i]);
            constantMatrix[i+2][3] = QualityUtils.qualToErrorProb(overallGCP[i]);
            constantMatrix[i+2][4] = QualityUtils.qualToErrorProb(deletionGOP[i]);
            constantMatrix[i+2][5] = QualityUtils.qualToErrorProb(overallGCP[i]);
        }
        constantMatrix[l+1][4] = 1.0;
        constantMatrix[l+1][5] = 1.0;

        // note that we initialized the constants
        constantsAreInitialized = true;
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
    private void updateCell( final int indI, final int indJ, final double prior, final double[] constants,
                             final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        matchMetricArray[indI][indJ] = prior * ( matchMetricArray[indI - 1][indJ - 1] * constants[0] +
                XMetricArray[indI - 1][indJ - 1] * constants[1] +
                YMetricArray[indI - 1][indJ - 1] * constants[1] );
        XMetricArray[indI][indJ] = matchMetricArray[indI - 1][indJ] * constants[2] + XMetricArray[indI - 1][indJ] * constants[3];
        YMetricArray[indI][indJ] = matchMetricArray[indI][indJ - 1] * constants[4] + YMetricArray[indI][indJ - 1] * constants[5];
    }
}
