using System;
using Bio;
using Bio.Math;

namespace Bio.PairHMM
{

	


	/// <summary>
	/// Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
	/// 
	/// User: rpoplin
	/// Date: 10/16/12
	/// </summary>
	public abstract class PairHMM
	{
		protected internal static readonly byte? MAX_CACHED_QUAL = byte.MaxValue;
		protected internal static readonly byte DEFAULT_GOP = (byte) 45;
		protected internal static readonly byte DEFAULT_GCP = (byte) 10;


		protected internal double[][] transition = null; // The transition probabilities cache
		protected internal double[][] prior = null; // The prior probabilities cache
		protected internal bool constantsAreInitialized = false;

		protected internal byte[] previousHaplotypeBases;



		public enum HMM_IMPLEMENTATION
		{
			/* Very slow implementation which uses very accurate log10 sum functions. Only meant to be used as a reference test implementation */
			EXACT,
			/* PairHMM as implemented for the UnifiedGenotyper. Uses log10 sum functions accurate to only 1E-4 */
			ORIGINAL,
			/* Optimized version of the PairHMM which caches per-read computations and operations in real space to avoid costly sums of log10'ed likelihoods */
			LOGLESS_CACHING
		}

		protected internal double[][] matchMatrix = null;
		protected internal double[][] insertionMatrix = null;
		protected internal double[][] deletionMatrix = null;
		protected internal int maxHaplotypeLength, maxReadLength;
		protected internal int paddedMaxReadLength, paddedMaxHaplotypeLength;
		protected internal int paddedReadLength, paddedHaplotypeLength;
		private bool initialized = false;

		/// <summary>
		/// Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
		/// 
		/// Note: Do not worry about padding, just provide the true max length of the read and haplotype. The HMM will take care of the padding.
		/// </summary>
		/// <param name="haplotypeMaxLength"> the max length of haplotypes we want to use with this PairHMM </param>
		/// <param name="readMaxLength"> the max length of reads we want to use with this PairHMM </param>
		public virtual void initialize(int readMaxLength, int haplotypeMaxLength)
		{
			if (readMaxLength <= 0)
			{
				throw new System.ArgumentException("READ_MAX_LENGTH must be > 0 but got " + readMaxLength);
			}
			if (haplotypeMaxLength <= 0)
			{
				throw new System.ArgumentException("HAPLOTYPE_MAX_LENGTH must be > 0 but got " + haplotypeMaxLength);
			}

			maxHaplotypeLength = haplotypeMaxLength;
			maxReadLength = readMaxLength;

			// M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
			paddedMaxReadLength = readMaxLength + 1;
			paddedMaxHaplotypeLength = haplotypeMaxLength + 1;


			matchMatrix = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, paddedMaxHaplotypeLength);
			insertionMatrix = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, paddedMaxHaplotypeLength);
			deletionMatrix = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, paddedMaxHaplotypeLength);

			previousHaplotypeBases = null;
			constantsAreInitialized = false;
			initialized = true;
		}



		/// <summary>
		/// Compute the total probability of read arising from haplotypeBases given base substitution, insertion, and deletion
		/// probabilities.
		/// 
		/// Note on using hapStartIndex.  This allows you to compute the exact true likelihood of a full haplotypes
		/// given a read, assuming that the previous calculation read over a full haplotype, recaching the read values,
		/// starting only at the place where the new haplotype bases and the previous haplotype bases different.  This
		/// index is 0-based, and can be computed with findFirstPositionWhereHaplotypesDiffer given the two haplotypes.
		/// Note that this assumes that the read and all associated quals values are the same.
		/// </summary>
		/// <param name="haplotypeBases"> the full sequence (in standard SAM encoding) of the haplotype, must be >= than read bases in length </param>
		/// <param name="readBases"> the bases (in standard encoding) of the read, must be <= haplotype bases in length </param>
		/// <param name="readQuals"> the phred-scaled per base substitition quality scores of read.  Must be the same length as readBases </param>
		/// <param name="insertionGOP"> the phred-scaled per base insertion quality scores of read.  Must be the same length as readBases </param>
		/// <param name="deletionGOP"> the phred-scaled per base deletion quality scores of read.  Must be the same length as readBases </param>
		/// <param name="overallGCP"> the phred-scaled gap continuation penalties scores of read.  Must be the same length as readBases </param>
		/// <param name="hapStartIndex"> start the hmm calculation at this offset in haplotype bases.  Used in the caching calculation
		///                      where multiple haplotypes are used, and they only diff starting at hapStartIndex </param>
		/// <param name="recacheReadValues"> if false, we don't recalculate any cached results, assuming that readBases and its associated
		///                          parameters are the same, and only the haplotype bases are changing underneath us </param>
		/// <returns> the log10 probability of read coming from the haplotype under the provided error model </returns>
		public double computeReadLikelihoodGivenHaplotypeLog10(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, int hapStartIndex, bool recacheReadValues)
		{
			if (!initialized)
			{
				throw new Exception("Must call initialize before calling computeReadLikelihoodGivenHaplotypeLog10");
			}
			if (haplotypeBases == null)
			{
				throw new System.ArgumentException("haplotypeBases cannot be null");
			}
			if (haplotypeBases.Length > maxHaplotypeLength)
			{
				throw new System.ArgumentException("Haplotype bases is too long, got " + haplotypeBases.Length + " but max is " + maxHaplotypeLength);
			}
			if (readBases == null)
			{
				throw new System.ArgumentException("readBases cannot be null");
			}
			if (readBases.Length > maxReadLength)
			{
				throw new System.ArgumentException("readBases is too long, got " + readBases.Length + " but max is " + maxReadLength);
			}
			if (readQuals.Length != readBases.Length)
			{
				throw new System.ArgumentException("Read bases and read quals aren't the same size: " + readBases.Length + " vs " + readQuals.Length);
			}
			if (insertionGOP.Length != readBases.Length)
			{
				throw new System.ArgumentException("Read bases and read insertion quals aren't the same size: " + readBases.Length + " vs " + insertionGOP.Length);
			}
			if (deletionGOP.Length != readBases.Length)
			{
				throw new System.ArgumentException("Read bases and read deletion quals aren't the same size: " + readBases.Length + " vs " + deletionGOP.Length);
			}
			if (overallGCP.Length != readBases.Length)
			{
				throw new System.ArgumentException("Read bases and overall GCP aren't the same size: " + readBases.Length + " vs " + overallGCP.Length);
			}

			paddedReadLength = readBases.Length + 1;
			paddedHaplotypeLength = haplotypeBases.Length + 1;

			// override the haplotype start index in case the previous haplotype had different length or this is a new read
			hapStartIndex = (previousHaplotypeBases == null || haplotypeBases.Length != previousHaplotypeBases.Length || recacheReadValues) ? 0 : hapStartIndex;


			double result = subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);

			if (!MathUtils.goodLog10Probability(result))
			{
				throw new Exception("PairHMM Log Probability cannot be greater than 0: " + string.Format("haplotype: {0}, read: {1}, result: {2:F}"));//, Arrays.ToString(haplotypeBases), Arrays.ToString(readBases), result));
			}

			// Warning: Careful if using the PairHMM in parallel! (this update has to be taken care of).
			// Warning: This assumes no downstream modification of the haplotype bases (saves us from copying the array). It is okay for the haplotype caller and the Unified Genotyper.
			previousHaplotypeBases = haplotypeBases;

			return result;
		}

		/// <summary>
		/// To be overloaded by subclasses to actually do calculation for #computeReadLikelihoodGivenHaplotypeLog10
		/// </summary>
		protected internal abstract double subComputeReadLikelihoodGivenHaplotypeLog10(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, int hapStartIndex, bool recacheReadValues);

		/// <summary>
		/// Print out the core hmm matrices for debugging
		/// </summary>
		protected internal virtual void dumpMatrices()
		{
			dumpMatrix("matchMetricArray", matchMatrix);
			dumpMatrix("insertionMatrix", insertionMatrix);
			dumpMatrix("deletionMatrix", deletionMatrix);
		}

		/// <summary>
		/// Print out in a human readable form the matrix for debugging </summary>
		/// <param name="name"> the name of this matrix </param>
		/// <param name="matrix"> the matrix of values </param>
		private void dumpMatrix(string name, double[][] matrix)
		{

			Console.Write("%s%n", name);
			for (int i = 0; i < matrix.Length; i++)
			{
				Console.Write("\t{0}[{1:D}]", name, i);
				for (int j = 0; j < matrix[i].Length; j++)
				{
					if (double.IsInfinity(matrix[i][j]))
					{
						Console.Write(" {0,15}", string.Format("{0:F}", matrix[i][j]));
					}
					else
					{

						Console.Write(" % 15.5e", matrix[i][j]);
					}
				}
				Console.WriteLine();
			}
		}

		/// <summary>
		/// Compute the first position at which two haplotypes differ
		/// 
		/// If the haplotypes are exact copies of each other, returns the min length of the two haplotypes.
		/// </summary>
		/// <param name="haplotype1"> the first haplotype1 </param>
		/// <param name="haplotype2"> the second haplotype1 </param>
		/// <returns> the index of the first position in haplotype1 and haplotype2 where the byte isn't the same </returns>
		public static int findFirstPositionWhereHaplotypesDiffer(byte[] haplotype1, byte[] haplotype2)
		{
			if (haplotype1 == null || haplotype1.Length == 0)
			{
				throw new System.ArgumentException("Haplotype1 is bad ");// + Arrays.ToString(haplotype1));
			}
			if (haplotype2 == null || haplotype2.Length == 0)
			{
				throw new System.ArgumentException("Haplotype2 is bad ");// + Arrays.ToString(haplotype2));
			}

			for (int iii = 0; iii < haplotype1.Length && iii < haplotype2.Length; iii++)
			{
				if (haplotype1[iii] != haplotype2[iii])
				{
					return iii;
				}
			}

			return System.Math.Min(haplotype1.Length, haplotype2.Length);
		}
	}

}