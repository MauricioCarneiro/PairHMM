package org.broadinstitute.pairhmm;

import org.broadinstitute.utils.QualityUtils;

/**
 * Features of this implementation over standard:
 *
 *  -- caching constants
 *  -- no logs throughout the routine (only the result)
 *  -- does not recalculate matrix for similar incoming haplotypes (uses hapStartIndex)
 *  -- uses single matrix for all operations
 *  -- uses path matrix to keep track of the state of the most likely path
 *  -- only calculates most likely path, does not sum over all possible paths
 *
 *  warning: this implementation has not yet been validated! Its experimental (read the name)
 *
 *
 * Created with IntelliJ IDEA.
 * User: rpoplin, carneiro
 * Date: 10/16/12
 */
public class ExperimentalPairHMM extends PairHMM {
    private double[][] transition = null; // The cache
    private double[][] prior = null; // The cache
    private int[][] path = null; // keeps track of the status of the path taken (to choose between match (0), unmatch (1), deletionGOP (2), deletionGCP(3), insertionGOP (4), insertionGCP (5))
    private boolean transitionMatrixInitialized = false;

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength) {
        // X and Y are size of read and haplotype + 1 because of an extra row and column for initial conditions (all 1's)
        X_METRIC_MAX_LENGTH = readMaxLength + 1;
        Y_METRIC_MAX_LENGTH = haplotypeMaxLength + 1;


        maxHaplotypeLength = haplotypeMaxLength;
        maxReadLength = readMaxLength;

        path = new int[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];                // correctly initializes everyone as a match
        transition = new double[X_METRIC_MAX_LENGTH][5];                         // will be initialized by every read
        prior = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];            // will be initialized by every haplotype and read (from the point where they first differ from the previous one)
        matchMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH]; // initialized here with a cusion row/column to give alignments starting in any position the same probability

        for (int j = 0; j < Y_METRIC_MAX_LENGTH; j++) {
            matchMetricArray[0][j] = 1.0;
        }

        initialized = true;
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
    public void initializePriorMatrix(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {

        // initialize the distance matrix for all combinations of read x haplotype bases
        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                final boolean isMatch = x == y || x == (byte) 'N' || y == (byte) 'N';
                prior[i+1][j+1] = isMatch ? QualityUtils.qualToProb(qual) : QualityUtils.qualToErrorProb(qual);
            }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * The matrix will hold an array for each base in the read with the following constants (given the index):
     *
     * 0 : probability of not being an insertion or a deletion (i.e. the probability of a match/mismatch)
     * 1 : probability of an insertion
     * 2 : probability of continuing inside an insertion
     * 3 : probability of a deletion
     * 4 : probability of continuing inside a deletion
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    private void initializeTransitionMatrix(final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        final int l = insertionGOP.length;

        for (int i = 0; i < l; i++) {
            final int qualIndexGOP = Math.min(insertionGOP[i] + deletionGOP[i], Byte.MAX_VALUE);
            transition[i+1][0] = QualityUtils.qualToProb((byte) qualIndexGOP);
            transition[i+1][1] = QualityUtils.qualToErrorProb(insertionGOP[i]);
            transition[i+1][2] = QualityUtils.qualToErrorProb(overallGCP[i]);
            transition[i+1][3] = QualityUtils.qualToErrorProb(deletionGOP[i]);
            transition[i+1][4] = QualityUtils.qualToErrorProb(overallGCP[i]);
        }
        transition[l][3] = 1.0; // allow deletions in the end of the read for free (note l is the last position in the read here since the matrix goes to l+1)
        transition[l][4] = 1.0; // allow deletion continuation in the end of the read for free

        transitionMatrixInitialized = true;
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
                                                               final boolean reinitializeTransitionMatrix ) {
        if (!transitionMatrixInitialized || reinitializeTransitionMatrix) {
            initializeTransitionMatrix(insertionGOP, deletionGOP, overallGCP);
        }
        initializePriorMatrix(haplotypeBases, readBases, readQuals, hapStartIndex);

        // we need to only operate over X and Y according to this read and haplotype lengths, not the max lengths
        final int readXMetricLength = readBases.length + 1;
        final int hapYMetricLength = haplotypeBases.length + 1;

        double maxLikelihood = -1.0;
        for (int i = 1; i < readXMetricLength; i++) {
//            dumpMatrices();
            for (int j = hapStartIndex + 1; j < hapYMetricLength; j++) {
                matchMetricArray[i][j] = updateCell(i, j);
                maxLikelihood = maxLikelihood > matchMetricArray[i][j] ? maxLikelihood : matchMetricArray[i][j];
            }
        }
        dumpMatrices();

        // final probability could be anywhere in the matrix
        return Math.log10(maxLikelihood);
    }


    /**
     * Updates a cell in the HMM matrix
     *
     *  match :    i-1 , j-1
     *  insertion: i-1 , j
     *  deletion:  i   , j-1
     *
     * @param i row index in the matrices to update
     * @param j column index in the matrices to update
     */
    private double updateCell(final int i, final int j) {
        final double match = prior[i][j] * matchMetricArray[i-1][j-1] * transition[i][0];           // prior * prob_getting_here * prob_of_(mis)match

        final int insertionPathIndex = (path[i-1][j] == 1 || path[i-1][j] == 2) ? 2 : 1;            // check if the path leading up to here was already an insertion, if it is, apply GCP for all other paths apply GOP.
        final double insertion = matchMetricArray[i-1][j] * transition[i][insertionPathIndex];      // prob_getting_here * prob_of_insertion/prob_of_ins_continuation

        final int deletionPathIndex = (path[i][j-1] == 3 || path[i][j-1] == 4) ? 4 : 3;             // check if the path leading up to here was already a deletion, if it is, apply GCP for all other paths apply GOP.
        final double deletion = matchMetricArray[i][j-1] * transition[i][deletionPathIndex];        // prob_getting_here * prob_of_deletion/prob_of_del_continuation

        double result;
        if (match >= insertion && match >= deletion) {
            path[i][j] = 0;
            result = match;
        } else if (insertion >= match && insertion >= deletion) {
            path[i][j] = insertionPathIndex;
            result = insertion;
        } else {
            path[i][j] = deletionPathIndex;
            result = deletion;
        }
        return result;
    }

    @Override
    protected void dumpMatrices() {
        dumpMatrix("matchMetricArray", matchMetricArray);
    }
}
