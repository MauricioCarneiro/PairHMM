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

package org.broadinstitute.utils;

/**
 * QualityUtils is a static class (no instantiation allowed!) with some utility methods for manipulating
 * quality scores.
 *
 * @author Kiran Garimella, Mark DePristo
 * @since Way back
 */
public class QualityUtils {

    /**
     * The lowest quality score for a base that is considered reasonable for statistical analysis.  This is
     * because Q 6 => you stand a 25% of being right, which means all bases are equally likely
     */
    public final static byte MIN_USABLE_Q_SCORE = 6;

    /**
     * Cached values for qual as byte calculations so they are very fast
     */
    private static float qualToErrorProbCache[] = new float[256];
    private static float qualToProbLog10Cache[] = new float[256];

    static {
        for (int i = 0; i < 256; i++) {
            qualToErrorProbCache[i] = qualToErrorProb((float) i);
            qualToProbLog10Cache[i] = (float) Math.log10(1.0 - qualToErrorProbCache[i]);
        }
    }

    /**
     * Private constructor.  No instantiating this class!
     */
    private QualityUtils() {}

    // ----------------------------------------------------------------------
    //
    // These are all functions to convert a phred-scaled quality score to a probability
    //
    // ----------------------------------------------------------------------

    /**
     * Convert a phred-scaled quality score to its probability of being true (Q30 => 0.999)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a discretized byte value, this function uses a cache so is very efficient
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a quality score (0-255)
     * @return a probability (0.0-1.0)
     */
    public static float qualToProb(final byte qual) {
        return 1.f - qualToErrorProb(qual);
    }

    /**
     * Convert a phred-scaled quality score to its log10 probability of being true (Q30 => log10(0.999))
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a float value, this function must call Math.pow so can be quite expensive
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a phred-scaled quality score encoded as a float.  Can be non-integer values (30.5)
     * @return a probability (0.0-1.0)
     */
    public static float qualToProbLog10(final byte qual) {
        return qualToProbLog10Cache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    /**
     * Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a float value, this function must call Math.pow so can be quite expensive
     *
     * @param qual a phred-scaled quality score encoded as a float.  Can be non-integer values (30.5)
     * @return a probability (0.0-1.0)
     */
    public static float qualToErrorProb(final float qual) {
        if ( qual < 0.0 ) throw new IllegalArgumentException("qual must be >= 0.0 but got " + qual);
        return (float) Math.pow(10.0, qual / -10.0);
    }

    /**
     * Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * Because the input is a byte value, this function uses a cache so is very efficient
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a phred-scaled quality score encoded as a byte
     * @return a probability (0.0-1.0)
     */
    public static float qualToErrorProb(final byte qual) {
        return qualToErrorProbCache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }


    /**
     * Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * The calculation is extremely efficient
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param qual a phred-scaled quality score encoded as a byte
     * @return a probability (0.0-1.0)
     */
    public static float qualToErrorProbLog10(final byte qual) {
        return qualToErrorProbLog10((float)(qual & 0xFF));
    }

    /**
     * Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * The calculation is extremely efficient
     *
     * @param qual a phred-scaled quality score encoded as a float
     * @return a probability (0.0-1.0)
     */
    public static float qualToErrorProbLog10(final float qual) {
        if ( qual < 0.0 ) throw new IllegalArgumentException("qual must be >= 0.0 but got " + qual);
        return qual / -10.f;
    }
}

