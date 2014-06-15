using System;
using Bio.Utils;
using Bio.Math;


namespace Bio.PairHMM
{

	

	public class LoglessPairHMM : PairHMM
	{
		protected internal static readonly double INITIAL_CONDITION = System.Math.Pow(2, 1020);
		protected internal static readonly double INITIAL_CONDITION_LOG10 = System.Math.Log10(INITIAL_CONDITION);

		/// <summary>
		/// {@inheritDoc}
		/// </summary>
		public override void initialize(int readMaxLength, int haplotypeMaxLength)
		{
			base.initialize(readMaxLength, haplotypeMaxLength);

			transition = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, 6);
            //ORIGINAL LINE: prior = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
			prior = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, paddedMaxHaplotypeLength);
		}

		/// <summary>
		/// {@inheritDoc}
		/// </summary>
		protected internal override double subComputeReadLikelihoodGivenHaplotypeLog10(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, int hapStartIndex, bool recacheReadValues)
		{

			if (previousHaplotypeBases == null || previousHaplotypeBases.Length != haplotypeBases.Length)
			{
				double initialValue = INITIAL_CONDITION / haplotypeBases.Length;
				// set the initial value (free deletions in the beginning) for the first row in the deletion matrix
				for (int j = 0; j < paddedHaplotypeLength; j++)
				{
					deletionMatrix[0][j] = initialValue;
				}
			}

			if (!constantsAreInitialized || recacheReadValues)
			{
				initializeProbabilities(insertionGOP, deletionGOP, overallGCP);
			}
			initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);

			for (int i = 1; i < paddedReadLength; i++)
			{
				// +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
				for (int j = hapStartIndex + 1; j < paddedHaplotypeLength; j++)
				{
					updateCell(i, j, prior[i][j], transition[i]);
                }
			}

			// final probability is the log10 sum of the last element in the Match and Insertion state arrays
			// this way we ignore all paths that ended in deletions! (huge)
			// but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
			int endI = paddedReadLength - 1;
			double finalSumProbabilities = 0.0;
			for (int j = 1; j < paddedHaplotypeLength; j++)
			{
				finalSumProbabilities += matchMatrix[endI][j] + insertionMatrix[endI][j];
			}
            return System.Math.Log10(finalSumProbabilities) - INITIAL_CONDITION_LOG10;
		}

		/// <summary>
		/// Initializes the matrix that holds all the constants related to the editing
		/// distance between the read and the haplotype.
		/// </summary>
		/// <param name="haplotypeBases"> the bases of the haplotype </param>
		/// <param name="readBases">      the bases of the read </param>
		/// <param name="readQuals">      the base quality scores of the read </param>
		/// <param name="startIndex">     where to start updating the distanceMatrix (in case this read is similar to the previous read) </param>
		public virtual void initializePriors(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, int startIndex)
		{

			// initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
			// Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

			for (int i = 0; i < readBases.Length; i++)
			{

				byte x = readBases[i];
				byte qual = readQuals[i];
				for (int j = startIndex; j < haplotypeBases.Length; j++)
				{
					byte y = haplotypeBases[j];
					prior[i + 1][j + 1] = (x == y || x == (byte) 'N' || y == (byte) 'N' ? QualityUtils.qualToProb(qual) : QualityUtils.qualToErrorProb(qual));
				}
			}
		}

		/// <summary>
		/// Initializes the matrix that holds all the constants related to quality scores.
		/// </summary>
		/// <param name="insertionGOP">   insertion quality scores of the read </param>
		/// <param name="deletionGOP">    deletion quality scores of the read </param>
		/// <param name="overallGCP">     overall gap continuation penalty </param>
		private void initializeProbabilities(byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP)
		{
			for (int i = 0; i < insertionGOP.Length; i++)
			{
				int qualIndexGOP = System.Math.Min(insertionGOP[i] + deletionGOP[i], byte.MaxValue);
				transition[i + 1][0] = QualityUtils.qualToProb((byte) qualIndexGOP);
				transition[i + 1][1] = QualityUtils.qualToProb(overallGCP[i]);
				transition[i + 1][2] = QualityUtils.qualToErrorProb(insertionGOP[i]);
				transition[i + 1][3] = QualityUtils.qualToErrorProb(overallGCP[i]);
				transition[i + 1][4] = QualityUtils.qualToErrorProb(deletionGOP[i]);
				transition[i + 1][5] = QualityUtils.qualToErrorProb(overallGCP[i]);
			}
			// note that we initialized the constants
			constantsAreInitialized = true;
		}

		/// <summary>
		/// Updates a cell in the HMM matrix
		/// 
		/// The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
		/// initial conditions
		/// </summary>
		/// <param name="indI">             row index in the matrices to update </param>
		/// <param name="indJ">             column index in the matrices to update </param>
		/// <param name="prior">            the likelihood editing distance matrix for the read x haplotype </param>
		/// <param name="transition">        an array with the six transition relevant to this location </param>
		private void updateCell(int indI, int indJ, double prior, double[] transition)
		{

			matchMatrix[indI][indJ] = prior * (matchMatrix[indI - 1][indJ - 1] * transition[0] + insertionMatrix[indI - 1][indJ - 1] * transition[1] + deletionMatrix[indI - 1][indJ - 1] * transition[1]);
			insertionMatrix[indI][indJ] = matchMatrix[indI - 1][indJ] * transition[2] + insertionMatrix[indI - 1][indJ] * transition[3];
			deletionMatrix[indI][indJ] = matchMatrix[indI][indJ - 1] * transition[4] + deletionMatrix[indI][indJ - 1] * transition[5];
		}
	}

}