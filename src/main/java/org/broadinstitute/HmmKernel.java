package org.broadinstitute;

import edu.syr.pcpratts.rootbeer.runtime.Kernel;
import org.broadinstitute.pairhmm.StandardPairHMM;
import org.broadinstitute.utils.QualityUtils;

public class HmmKernel implements Kernel {

    private double likelihood;

    private byte[] haplotypeBases;
    private byte[] readBases;
    private byte[] readQuals;
    private byte[] insertionGOP;
    private byte[] deletionGOP;
    private byte[] overallGCP;
    private int hapStartIndex;
    private boolean recacheReadValues;

    private static final int X_METRIC_LENGTH = 10000;
    private static final int Y_METRIC_LENGTH = 10000;

    public HmmKernel(final byte[] haplotypeBases, final byte[] readBases, 
      final byte[] readQuals, final byte[] insertionGOP, 
      final byte[] deletionGOP, final byte[] overallGCP, 
      final int hapStartIndex, final boolean recacheReadValues) {

        this.haplotypeBases = haplotypeBases;
        this.readBases = readBases;
        this.readQuals = readQuals;
        this.insertionGOP = insertionGOP;
        this.deletionGOP = deletionGOP;
        this.overallGCP = overallGCP;
        this.hapStartIndex = hapStartIndex;
        this.recacheReadValues = recacheReadValues; 
        this.likelihood = likelihood;
    }

    public void gpuMethod(){
      StandardPairHMM hmm = new StandardPairHMM();
      hmm.initialize(X_METRIC_LENGTH + 2, Y_METRIC_LENGTH + 2);
      likelihood = hmm.subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, 
        readBases, cleanupQualityScores(readQuals), insertionGOP, deletionGOP, 
        overallGCP, hapStartIndex, recacheReadValues);
    }
  
    public double getLikelihood(){
      return likelihood;
    }


    /**
     * Ensures that all the qual scores have valid values
     *
     * @param readQuals read qualities array (modified in place)
     * @return the readQuals array for convenience
     */
    private byte[] cleanupQualityScores(byte[] readQuals) {
        for (int i = 0; i < readQuals.length; i++) {
            readQuals[i] = (readQuals[i] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[i] > Byte.MAX_VALUE ? Byte.MAX_VALUE : readQuals[i]));
        }
        return readQuals;
    }
}
