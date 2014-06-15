using System;



namespace Bio.Utils
{

	/// <summary>
	/// QualityUtils is a static class (no instantiation allowed!) with some utility methods for manipulating
	/// quality scores.
	/// 
	/// @author Kiran Garimella, Mark DePristo
	/// @since Way back
	/// </summary>
	public class QualityUtils
	{

		/// <summary>
		/// The lowest quality score for a base that is considered reasonable for statistical analysis.  This is
		/// because Q 6 => you stand a 25% of being right, which means all bases are equally likely
		/// </summary>
		public const byte MIN_USABLE_Q_SCORE = 6;

		/// <summary>
		/// Cached values for qual as byte calculations so they are very fast
		/// </summary>
		private static double[] qualToErrorProbCache = new double[256];
		private static double[] qualToProbLog10Cache = new double[256];

		static QualityUtils()
		{
			for (int i = 0; i < 256; i++)
			{
				qualToErrorProbCache[i] = qualToErrorProb((double) i);
                qualToProbLog10Cache[i] = System.Math.Log10(1.0 - qualToErrorProbCache[i]);
			}
		}

		/// <summary>
		/// Private constructor.  No instantiating this class!
		/// </summary>
		private QualityUtils()
		{
		}

		// ----------------------------------------------------------------------
		//
		// These are all functions to convert a phred-scaled quality score to a probability
		//
		// ----------------------------------------------------------------------

		/// <summary>
		/// Convert a phred-scaled quality score to its probability of being true (Q30 => 0.999)
		/// 
		/// This is the Phred-style conversion, *not* the Illumina-style conversion.
		/// 
		/// Because the input is a discretized byte value, this function uses a cache so is very efficient
		/// 
		/// WARNING -- because this function takes a byte for maxQual, you must be careful in converting
		/// integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
		/// </summary>
		/// <param name="qual"> a quality score (0-255) </param>
		/// <returns> a probability (0.0-1.0) </returns>
		public static double qualToProb(byte qual)
		{
			return 1.0 - qualToErrorProb(qual);
		}

		/// <summary>
		/// Convert a phred-scaled quality score to its log10 probability of being true (Q30 => log10(0.999))
		/// 
		/// This is the Phred-style conversion, *not* the Illumina-style conversion.
		/// 
		/// Because the input is a double value, this function must call Math.pow so can be quite expensive
		/// 
		/// WARNING -- because this function takes a byte for maxQual, you must be careful in converting
		/// integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
		/// </summary>
		/// <param name="qual"> a phred-scaled quality score encoded as a double.  Can be non-integer values (30.5) </param>
		/// <returns> a probability (0.0-1.0) </returns>
		public static double qualToProbLog10(byte qual)
		{
			return qualToProbLog10Cache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
		}

		/// <summary>
		/// Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
		/// 
		/// This is the Phred-style conversion, *not* the Illumina-style conversion.
		/// 
		/// Because the input is a double value, this function must call Math.pow so can be quite expensive
		/// </summary>
		/// <param name="qual"> a phred-scaled quality score encoded as a double.  Can be non-integer values (30.5) </param>
		/// <returns> a probability (0.0-1.0) </returns>
		public static double qualToErrorProb(double qual)
		{
			if (qual < 0.0)
			{
				throw new System.ArgumentException("qual must be >= 0.0 but got " + qual);
			}
            return System.Math.Pow(10.0, qual / -10.0);
		}

		/// <summary>
		/// Convert a phred-scaled quality score to its probability of being wrong (Q30 => 0.001)
		/// 
		/// This is the Phred-style conversion, *not* the Illumina-style conversion.
		/// 
		/// Because the input is a byte value, this function uses a cache so is very efficient
		/// 
		/// WARNING -- because this function takes a byte for maxQual, you must be careful in converting
		/// integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
		/// </summary>
		/// <param name="qual"> a phred-scaled quality score encoded as a byte </param>
		/// <returns> a probability (0.0-1.0) </returns>
		public static double qualToErrorProb(byte qual)
		{
			return qualToErrorProbCache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
		}


		/// <summary>
		/// Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
		/// 
		/// This is the Phred-style conversion, *not* the Illumina-style conversion.
		/// 
		/// The calculation is extremely efficient
		/// 
		/// WARNING -- because this function takes a byte for maxQual, you must be careful in converting
		/// integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
		/// </summary>
		/// <param name="qual"> a phred-scaled quality score encoded as a byte </param>
		/// <returns> a probability (0.0-1.0) </returns>
		public static double qualToErrorProbLog10(byte qual)
		{
			return qualToErrorProbLog10((double)(qual & 0xFF));
		}

		/// <summary>
		/// Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
		/// 
		/// This is the Phred-style conversion, *not* the Illumina-style conversion.
		/// 
		/// The calculation is extremely efficient
		/// </summary>
		/// <param name="qual"> a phred-scaled quality score encoded as a double </param>
		/// <returns> a probability (0.0-1.0) </returns>
		public static double qualToErrorProbLog10(double qual)
		{
			if (qual < 0.0)
			{
				throw new System.ArgumentException("qual must be >= 0.0 but got " + qual);
			}
			return qual / -10.0;
		}
	}


}