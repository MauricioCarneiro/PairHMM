package org.broadinstitute;

import edu.syr.pcpratts.rootbeer.runtime.Kernel;

public class HmmKernel implements Kernel {

    private double result;
    private double likelihood;

    private byte[] haplotypeBases;
    private byte[] readBases;
    private byte[] readQuals;
    private byte[] insertionGOP;
    private byte[] deletionGOP;
    private byte[] overallGCP;
    private int hapStartIndex;
    private boolean recacheReadValues;

    public HmmKernel(final byte[] haplotypeBases, final byte[] readBases, 
      final byte[] readQuals, final byte[] insertionGOP, 
      final byte[] deletionGOP, final byte[] overallGCP, 
      final int hapStartIndex, final boolean recacheReadValues,
      final double likelihood) {

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
      result = likelihood;
    }

    public double getResult(){
      return result;
    }
                        
    public double getLikelihood(){
      return likelihood;
    }
}
