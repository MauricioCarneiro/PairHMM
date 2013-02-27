package org.broadinstitute;

/**
 * Created with IntelliJ IDEA.
 * User: tjordan
 * Date: 1/14/13
 * Time: 2:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class TestRow {
    private byte[] haplotypeBases;
    private byte[] readBases;
    private byte[] readQuals;
    private byte[] readInsQuals;
    private byte[] readDelQuals;
    private byte[] overallGCP;
    private int haplotypeStart;
    private boolean reachedReadValue;
    private Double likelihood;

    TestRow(byte[] haplotypeBases, byte[] readBases,
            byte[] readQuals, byte[] readInsQuals,
            byte[] readDelQuals, byte[] overallGCP,
            int haplotypeStart, boolean reachedReadValue,
            Double likelihood) {

        this.haplotypeBases = haplotypeBases;
        this.readBases = readBases;
        this.readQuals = readQuals;
        this.readInsQuals = readInsQuals;
        this.readDelQuals = readDelQuals;
        this.overallGCP = overallGCP;
        this.haplotypeStart = haplotypeStart;
        this.reachedReadValue = reachedReadValue;
        this.likelihood = likelihood;
    }

    public byte[] getHaplotypeBases() {
        return this.haplotypeBases;
    }

    public byte[] getReadBases() {
        return this.readBases;
    }

    public byte[] getReadQuals() {
        return this.readQuals;
    }

    public byte[] getReadInsQuals() {
        return this.readInsQuals;
    }

    public byte[] getReadDelQuals() {
        return this.readDelQuals;
    }

    public byte[] getOverallGCP() {
        return this.overallGCP;
    }

    public int getHaplotypeStart() {
        return this.haplotypeStart;
    }

    public boolean getReachedReadValye() {
        return this.reachedReadValue;
    }

    public double getLikelihood() {
        return this.likelihood;
    }

}
