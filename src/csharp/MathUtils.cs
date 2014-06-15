using System;

namespace Bio.Math
{

	/// <summary>
	/// @author Mauricio Carneiro
	/// @since 8/16/12
	/// </summary>
	public class MathUtils
	{

		private static readonly double[] jacobianLogTable;
		private const double JACOBIAN_LOG_TABLE_STEP = 0.001;
		private const double MAX_JACOBIAN_TOLERANCE = 8.0;
		private const double JACOBIAN_LOG_TABLE_INV_STEP = 1.0 / 0.001;
		private static readonly int JACOBIAN_LOG_TABLE_SIZE = (int)(MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1;

		static MathUtils()
		{
			jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];

			for (int k = 0; k < JACOBIAN_LOG_TABLE_SIZE; k++)
			{
                jacobianLogTable[k] = System.Math.Log10(1.0 + System.Math.Pow(10.0, -((double)k) * JACOBIAN_LOG_TABLE_STEP));

			}
		}

		/// <summary>
		/// Checks that the result is a well-formed log10 probability
		/// </summary>
		/// <param name="result"> a supposedly well-formed log10 probability value.  By default allows
		///               -Infinity values, as log10(0.0) == -Infinity. </param>
		/// <returns> true if result is really well formed </returns>
		public static bool goodLog10Probability(double result)
		{
			return goodLog10Probability(result, true);
		}

		/// <summary>
		/// Checks that the result is a well-formed log10 probability
		/// </summary>
		/// <param name="result"> a supposedly well-formed log10 probability value </param>
		/// <param name="allowNegativeInfinity"> should we consider a -Infinity value ok? </param>
		/// <returns> true if result is really well formed </returns>
		public static bool goodLog10Probability(double result, bool allowNegativeInfinity)
		{
			return result <= 0.0 && result != double.PositiveInfinity && (allowNegativeInfinity || result != double.NegativeInfinity) && !double.IsNaN(result);
		}

		// A fast implementation of the Math.round() method.  This method does not perform
		// under/overflow checking, so this shouldn't be used in the general case (but is fine
		// if one is already make those checks before calling in to the rounding).
		public static int fastRound(double d)
		{
			return (d > 0.0) ? (int)(d + 0.5d) : (int)(d - 0.5d);
		}

		public static int maxElementIndex(double[] array, int endIndex)
		{
			if (array == null || array.Length == 0)
			{
				throw new System.ArgumentException("Array cannot be null!");
			}

			int maxI = 0;
			for (int i = 1; i < endIndex; i++)
			{
				if (array[i] > array[maxI])
				{
					maxI = i;
				}
			}

			return maxI;
		}

		public static int maxElementIndex(double[] array)
		{
			return maxElementIndex(array, array.Length);
		}

		public static double arrayMax(double[] array, int endIndex)
		{
			return array[maxElementIndex(array, endIndex)];
		}

		public static double log10sumLog10(double[] log10p, int start, int finish)
		{
			double sum = 0.0;

			double maxValue = arrayMax(log10p, finish);
			if (maxValue == double.NegativeInfinity)
			{
				return maxValue;
			}

			for (int i = start; i < finish; i++)
			{
                sum += System.Math.Pow(10.0, log10p[i] - maxValue);
			}

            return System.Math.Log10(sum) + maxValue;
		}

		public static double approximateLog10SumLog10(double[] vals)
		{
			return approximateLog10SumLog10(vals, vals.Length);
		}

		public static double approximateLog10SumLog10(double[] vals, int endIndex)
		{

//JAVA TO C# CONVERTER WARNING: The original Java variable was marked 'final':
//ORIGINAL LINE: final int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
			int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
			double approxSum = vals[maxElementIndex];

			for (int i = 0; i < endIndex; i++)
			{
				if (i == maxElementIndex || vals[i] == double.NegativeInfinity)
				{
					continue;
				}

//JAVA TO C# CONVERTER WARNING: The original Java variable was marked 'final':
//ORIGINAL LINE: final double diff = approxSum - vals[i];
				double diff = approxSum - vals[i];
				if (diff < MathUtils.MAX_JACOBIAN_TOLERANCE)
				{
					// See notes from the 2-inout implementation below
//JAVA TO C# CONVERTER WARNING: The original Java variable was marked 'final':
//ORIGINAL LINE: final int ind = fastRound(diff * MathUtils.JACOBIAN_LOG_TABLE_INV_STEP);
					int ind = fastRound(diff * MathUtils.JACOBIAN_LOG_TABLE_INV_STEP); // hard rounding
					approxSum += MathUtils.jacobianLogTable[ind];
				}
			}

			return approxSum;
		}

		public static double log10sumLog10(double[] log10p, int start)
		{
			return log10sumLog10(log10p, start, log10p.Length);
		}

		public static double log10sumLog10(double[] log10values)
		{
			return log10sumLog10(log10values, 0);
		}

		public static double approximateLog10SumLog10(double small, double big)
		{
			// make sure small is really the smaller value
			if (small > big)
			{
//JAVA TO C# CONVERTER WARNING: The original Java variable was marked 'final':
//ORIGINAL LINE: final double t = big;
				double t = big;
				big = small;
				small = t;
			}

			if (small == double.NegativeInfinity)
			{
				return big;
			}

//JAVA TO C# CONVERTER WARNING: The original Java variable was marked 'final':
//ORIGINAL LINE: final double diff = big - small;
			double diff = big - small;
			if (diff >= MAX_JACOBIAN_TOLERANCE)
			{
				return big;
			}

			// OK, so |y-x| < tol: we use the following identity then:
			// we need to compute log10(10^x + 10^y)
			// By Jacobian logarithm identity, this is equal to
			// max(x,y) + log10(1+10^-abs(x-y))
			// we compute the second term as a table lookup with integer quantization
			// we have pre-stored correction for 0,0.1,0.2,... 10.0
//JAVA TO C# CONVERTER WARNING: The original Java variable was marked 'final':
//ORIGINAL LINE: final int ind = fastRound(diff * JACOBIAN_LOG_TABLE_INV_STEP);
			int ind = fastRound(diff * JACOBIAN_LOG_TABLE_INV_STEP); // hard rounding
			return big + jacobianLogTable[ind];
		}
	}

    public static class Arrays
    {
        public static void fill(double[] arr, double val)
        {
            for(int i=0;i<arr.Length;i++)
            {
                arr[i]=val;
            }
        }
    }


}