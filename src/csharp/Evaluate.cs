using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Linq;

namespace org.broadinstitute
{

	using Log10PairHMM = Bio.PairHMM.Log10PairHMM;
    using LoglessPairHMM = Bio.PairHMM.LoglessPairHMM;
    using PairHMM = Bio.PairHMM.PairHMM;
    using QualityUtils = Bio.Utils.QualityUtils;

	/// <summary>
	/// Evaluating routine for the PairHMM
	/// 
	/// @author Mauricio Carneiro
	/// @since 8/15/12
	/// </summary>
    public class Evaluate
    {
        public readonly IList<PairHMM> pairHMM = new List<PairHMM>(20);
        private const int X_METRIC_LENGTH = 400;
        private const int Y_METRIC_LENGTH = 400;
    

        /// <summary>
        /// Initializes and computes the Pair HMM matrix.
        /// <p/>
        /// Use this method if you're calculating the entire matrix from scratch.
        /// </summary>
        /// <param name="haplotypeBases"> reference sequence bases </param>
        /// <param name="readBases">      comparison haplotype bases </param>
        /// <param name="readQuals">      comparison haplotype base quals (phred-scaled) </param>
        /// <param name="insertionGOP">   comparison haplotype insertion quals (phred-scaled) </param>
        /// <param name="deletionGOP">    comparison haplotype deletion quals (phred-scaled) </param>
        /// <param name="overallGCP">     comparison haplotype gap continuation quals (phred-scaled) </param>
        /// <returns> the likelihood of the alignment between read and haplotype </returns>
        public static double runhmm(PairHMM hmm, byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] insertionGOP, byte[] deletionGOP, byte[] overallGCP, int hapStartIndex, bool recacheReadValues)
        {
            return hmm.computeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, cleanupQualityScores(readQuals), insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);
        }


        private static string createFileName(string hmmName)
        {
            int runNumber = 0;
            string filename;

            do
            {
                runNumber++;
                filename = hmmName + "." + runNumber + ".out";
            } while (File.Exists(filename));
            return filename;
        }

        /// <summary>
        /// Ensures that all the qual scores have valid values
        /// </summary>
        /// <param name="readQuals"> read qualities array (modified in place) </param>
        /// <returns> the readQuals array for convenience </returns>
        private static byte[] cleanupQualityScores(byte[] readQuals)
        {
            for (int i = 0; i < readQuals.Length; i++)
            {
                readQuals[i] = (readQuals[i] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[i] > byte.MaxValue ? byte.MaxValue : readQuals[i]));
            }
            return readQuals;
        }

        private static void runTests(PairHMM hmm, IEnumerable<TestRow> testCache, string testSet)
        {
           
            string runName = "Test." + testSet;
            string filename = createFileName(runName);
            StreamWriter fout = new StreamWriter(filename);
            hmm.initialize(X_METRIC_LENGTH, Y_METRIC_LENGTH);
            foreach(var currentTest in testCache)           
            {            
                double likelihood = Evaluate.runhmm(hmm, currentTest.haplotypeBases, currentTest.readBases, currentTest.readQuals, currentTest.readInsQuals, currentTest.readDelQuals, currentTest.overallGCP, currentTest.haplotypeStart, currentTest.reachedReadValue);
                fout.Write("" + likelihood + "\n");
            }
            fout.Close();
        }


        private static IEnumerable<TestRow> createIteratorFor(string filePath)
        {
            return new ParserIterator(filePath).GetTestLines();
        }

        public static void Main(string[] argv)
        {

            var args = new HashSet<string>(argv);
            if (args.Count < 1)
            {
                throw new Exception("\r\nYou must specify a file name for input.\n" + "filename \n ----------------------------\n" + "Run with -d or --debug for debug output");
            }
            else
            {
                List<string> testFiles = new List<string>();
                foreach (string arg in args)
                {
                    if (!arg.StartsWith("-"))
                    {
                        testFiles.Add(arg);
                    }
                }

                foreach (String testSet in testFiles)
                {
                    var hmm = new Log10PairHMM(true);
                    var name = new FileInfo(testSet);
                    var testName = name.Name.Replace(name.Extension,"");
                    Evaluate.runTests(hmm, createIteratorFor(testSet), testName );
                    
                }
            }

        }

        /// <summary>
        /// Classes and Functions used to parse and apply the tests to the
        /// PairHMM class
        /// </summary>

        internal class ParserIterator 
        {

            IEnumerable<string> filelines;
            //private XReadLines m_reader;
            private TestRow m_nextEntry;
            public ParserIterator(string filePath)
            {
                //m_reader = new XReadLines(new File(filePath));
                filelines = File.ReadLines(filePath);
            }
            public IEnumerable<TestRow> GetTestLines()
            {
                foreach (var line in filelines)
                {
                    string[] p = line.Split(' ');
                    bool reachedRead = p[7].Trim().Equals("true");
                    yield return new TestRow(StringToBytes(p[0]), StringToBytes(p[1]), convertToByteArray(p[2]), convertToByteArray(p[3]), convertToByteArray(p[4]), convertToByteArray(p[5]), Convert.ToInt32(p[6]), reachedRead);

                    //      logger.debug("Haplotype Bases:" + p[0]);
                    //      logger.debug("Read Bases:" + p[1]);
                    //     logger.debug("Read Quals: " + p[2]);
                    //     logger.debug("Read Ins Quals: " + p[3]);
                    //     logger.debug("Read Del Quals: " + p[4]);
                    //     logger.debug("Overall GCP:" + p[5]);
                    //     logger.debug("Haplotype start:" + p[6]);
                    //      logger.debug("Boolean (reached read values):" + p[7]);

                }
            }
            private byte[] StringToBytes(string var)
            {
                return Encoding.ASCII.GetBytes(var);//.Select(x=>(byte)).to;
            }
#if FALSE
			public override bool hasNext()
			{
				if (!m_reader.hasNext())
				{
					return false;
				}

				string line = m_reader.next();
				string[] p = line.Split(" ", true);

				bool reachedRead = p[7].Trim().Equals("true");
				m_nextEntry = new TestRow(StringToBytes(p[0]), StringToBytes(p[1]), convertToByteArray(p[2]), convertToByteArray(p[3]), convertToByteArray(p[4]), convertToByteArray(p[5]), Convert.ToInt32(p[6]), reachedRead);

		  //      logger.debug("Haplotype Bases:" + p[0]);
		  //      logger.debug("Read Bases:" + p[1]);
		   //     logger.debug("Read Quals: " + p[2]);
		   //     logger.debug("Read Ins Quals: " + p[3]);
		   //     logger.debug("Read Del Quals: " + p[4]);
		   //     logger.debug("Overall GCP:" + p[5]);
		   //     logger.debug("Haplotype start:" + p[6]);
		  //      logger.debug("Boolean (reached read values):" + p[7]);

				return true;
			}
            
            public override TestRow next()
            {
                return m_nextEntry;
            }

            public override void remove()
            {
                throw new System.NotSupportedException("Not supported yet.");
            }
#endif
            private byte[] convertToByteArray(string quals)
            {
                byte[] output = StringToBytes(quals);

                for (int i = 0; i < output.Length; i++)
                {
                    output[i] -= 33;
                }

                return output;
            }

        }

        internal class TestRow
        {
            //Converted all bytes to byte
            internal byte[] haplotypeBases;
            internal byte[] readBases;
            internal byte[] readQuals;
            internal byte[] readInsQuals;
            internal byte[] readDelQuals;
            internal byte[] overallGCP;
            internal int haplotypeStart;
            internal bool reachedReadValue;

            internal TestRow(byte[] haplotypeBases, byte[] readBases, byte[] readQuals, byte[] readInsQuals, byte[] readDelQuals, byte[] overallGCP, int haplotypeStart, bool reachedReadValue)
            {
                this.haplotypeBases = haplotypeBases;
                this.readBases = readBases;
                this.readQuals = readQuals;
                this.readInsQuals = readInsQuals;
                this.readDelQuals = readDelQuals;
                this.overallGCP = overallGCP;
                this.haplotypeStart = haplotypeStart;
                this.reachedReadValue = reachedReadValue;
            }
        }
#if FALSE
		internal class XReadLines : IEnumerator<string>, IEnumerable<string>
		{
			private readonly BufferedReader @in; // The stream we're reading from
			private string nextLine = null; // Return value of next call to next()
			private readonly bool trimWhitespace;
			private readonly string commentPrefix;


			public XReadLines(File filename) : this(new FileReader(filename), true, null)
			{
			}

			/// <summary>
			/// Creates a new xReadLines object to read lines from an bufferedReader
			/// </summary>
			/// <param name="reader">         file name </param>
			/// <param name="trimWhitespace"> trim whitespace </param>
			/// <param name="commentPrefix">  prefix for comments or null if no prefix is set </param>
			public XReadLines(Reader reader, bool trimWhitespace, string commentPrefix)
			{
				this.@in = (reader is BufferedReader) ? (BufferedReader) reader : new BufferedReader(reader);
				this.trimWhitespace = trimWhitespace;
				this.commentPrefix = commentPrefix;
				try
				{
					this.nextLine = readNextLine();
				}
				catch (IOException e)
				{
					throw new System.ArgumentException(e);
				}
			}

			/// <summary>
			/// I'm an iterator too...
			/// </summary>
			/// <returns> an iterator </returns>
			public virtual IEnumerator<string> GetEnumerator()
			{
				return this;
			}

			public virtual bool hasNext()
			{
				return this.nextLine != null;
			}

			/// <summary>
			/// Actually reads the next line from the stream, not accessible publicly
			/// </summary>
			/// <returns> the next line or null </returns>
			/// <exception cref="java.io.IOException"> if an error occurs </exception>

			private string readNextLine()
			{
				string nextLine;
				while ((nextLine = this.@in.readLine()) != null)
				{
					if (this.trimWhitespace)
					{
						nextLine = nextLine.Trim();
						if (nextLine.Length == 0)
						{
							continue;
						}
					}
					if (this.commentPrefix != null)
					{
						if (nextLine.StartsWith(this.commentPrefix))
						{
							continue;
						}
					}
					break;
				}
				return nextLine;
			}

			/// <summary>
			/// Returns the next line (optionally minus whitespace)
			/// </summary>
			/// <returns> the next line </returns>
			public virtual string next()
			{
				try
				{
					string result = this.nextLine;
					this.nextLine = readNextLine();

					// If we haven't reached EOF yet
					if (this.nextLine == null)
					{
						@in.close(); // And close on EOF
					}

					// Return the line we read last time through.
					return result;
				}
				catch (IOException e)
				{
					throw new System.ArgumentException(e);
				}
			}

			// The file is read-only; we don't allow lines to be removed.
			public virtual void remove()
			{
				throw new System.NotSupportedException();
			}

		}
	}
#endif
    }
}