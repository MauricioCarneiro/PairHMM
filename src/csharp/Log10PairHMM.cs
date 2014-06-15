using System;
using Bio;
using Bio.Utils;
using Bio.Math;


namespace Bio.PairHMM
{

	/// <summary>
	/// Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
	/// </summary>
	public sealed class Log10PairHMM : PairHMM
	{		
		/// <summary>
		/// Create an uninitialized PairHMM
		/// </summary>
		/// <param name="doExactLog10"> should the log10 calculations be exact (slow) or approximate (faster) </param>
		public Log10PairHMM(bool doExactLog10)
		{
            DoingExactLog10Calculations = doExactLog10;
		}

		/// <summary>
		/// Is this HMM using exact log10 calculations? </summary>
		/// <returns> true if exact, false if approximate </returns>
		public bool DoingExactLog10Calculations
		{
			get; private set;
		}

		/// <summary>
		/// {@inheritDoc}
		/// </summary>
		public override void initialize(int readMaxLength, int haplotypeMaxLength)
		{
			base.initialize(readMaxLength, haplotypeMaxLength);

			for (int iii = 0; iii < paddedMaxReadLength; iii++)
			{
				Arrays.fill(matchMatrix[iii], double.NegativeInfinity);
				Arrays.fill(insertionMatrix[iii], double.NegativeInfinity);
				Arrays.fill(deletionMatrix[iii], double.NegativeInfinity);
			}

			transition = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, 6);
			prior = RectangularArrays.ReturnRectangularDoubleArray(paddedMaxReadLength, paddedMaxHaplotypeLength);
		}

        //The override method that calls this function, note that this probably should not be an override!!!!	
		protected internal override double subComputeReadLikelihoodGivenHaplotypeLog10(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, int hapStartIndex, bool recacheReadValues)
		{

			if (previousHaplotypeBases == null || previousHaplotypeBases.Length != haplotypeBases.Length)
			{
				// set the initial value (free deletions in the beginning) for the first row in the deletion matrix

				double initialValue = System.Math.Log10(1.0 / haplotypeBases.Length);
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
			double finalSumProbabilities = myLog10SumLog10(new double[]{matchMatrix[endI][1], insertionMatrix[endI][1]});
			for (int j = 2; j < paddedHaplotypeLength; j++)
			{
				finalSumProbabilities = myLog10SumLog10(new double[]{finalSumProbabilities, matchMatrix[endI][j], insertionMatrix[endI][j]});
			}

			return finalSumProbabilities;
		}

		/// <summary>
		/// Initializes the matrix that holds all the constants related to the editing
		/// distance between the read and the haplotype.
		/// </summary>
		/// <param name="haplotypeBases"> the bases of the haplotype </param>
		/// <param name="readBases">      the bases of the read </param>
		/// <param name="readQuals">      the base quality scores of the read </param>
		/// <param name="startIndex">     where to start updating the distanceMatrix (in case this read is similar to the previous read) </param>
		public void initializePriors(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, int startIndex)
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
					prior[i + 1][j + 1] = (x == y || x == (byte) 'N' || y == (byte) 'N' ? QualityUtils.qualToProbLog10(qual) : QualityUtils.qualToErrorProbLog10(qual));
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
				transition[i + 1][0] = QualityUtils.qualToProbLog10((byte) qualIndexGOP);
				transition[i + 1][1] = QualityUtils.qualToProbLog10(overallGCP[i]);
				transition[i + 1][2] = QualityUtils.qualToErrorProbLog10(insertionGOP[i]);
				transition[i + 1][3] = QualityUtils.qualToErrorProbLog10(overallGCP[i]);
				transition[i + 1][4] = QualityUtils.qualToErrorProbLog10(deletionGOP[i]);
				transition[i + 1][5] = QualityUtils.qualToErrorProbLog10(overallGCP[i]);
			}

			// note that we initialized the constants
			constantsAreInitialized = true;
		}


		/// <summary>
		/// Compute the log10SumLog10 of the values
		/// 
		/// NOTE NOTE NOTE
		/// 
		/// Log10PairHMM depends critically on this function tolerating values that are all -Infinity
		/// and the sum returning -Infinity.  Note good.  Needs to be fixed.
		/// 
		/// NOTE NOTE NOTE
		/// </summary>
		/// <param name="values"> an array of log10 probabilities that need to be summed </param>
		/// <returns> the log10 of the sum of the probabilities </returns>
		private double myLog10SumLog10(double[] values)
		{
			return DoingExactLog10Calculations ? MathUtils.log10sumLog10(values) : MathUtils.approximateLog10SumLog10(values);
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

			matchMatrix[indI][indJ] = prior + myLog10SumLog10(new double[]{matchMatrix[indI - 1][indJ - 1] + transition[0], insertionMatrix[indI - 1][indJ - 1] + transition[1], deletionMatrix[indI - 1][indJ - 1] + transition[1]});
			insertionMatrix[indI][indJ] = myLog10SumLog10(new double[] {matchMatrix[indI - 1][indJ] + transition[2], insertionMatrix[indI - 1][indJ] + transition[3]});
			deletionMatrix[indI][indJ] = myLog10SumLog10(new double[] {matchMatrix[indI][indJ - 1] + transition[4], deletionMatrix[indI][indJ - 1] + transition[5]});
		}
	}

}