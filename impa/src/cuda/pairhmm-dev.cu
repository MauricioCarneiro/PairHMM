#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>

//#include "cuPrintf.cu"

using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;

#define NOT_DONE 0
#define DONE 1
#define GAP_TO(sz, n) ((n - ((sz) % n)) % n)
#define PH2PR(n) (pow(NUM(10.0), -NUM(n) / NUM(10.0)))

#define MAX_HAPLEN 3072
#define NBLOCKS 150
#define NUM_OF_TESTCASES_PER_ITERATION 10000

template<class T> __device__ static inline T INITIAL_CONSTANT();
template<> __device__ static inline float INITIAL_CONSTANT<float>() { return 1e32f; }
template<> __device__ static inline double INITIAL_CONSTANT<double>() { return ldexp(1.0, 1020); }
template<class T> __device__ static inline T MIN_ACCEPTED();
template<> __device__ static inline float MIN_ACCEPTED<float>() { return 1e-28f; }
template<> __device__ static inline double MIN_ACCEPTED<double>() { return 0.0; }

struct Pair 
{
	double result;
	int offset_hap, offset_rs, offset_qidc, haplen, rslen, status;
};

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		cerr << "Error code: " << code << "\n";
		cerr << "Description: " << cudaGetErrorString(code) << "\n";
		cerr << "Line: " << line << "\n";
		exit(code);
	}
}

template <int nThreads, typename NUM> 
__global__ void compute_full_scores(char *g_chunk, int num_of_pairs, Pair *g_pair, int *g_pair_index, NUM *g_last_lines)
{

	// *******************************************************************************************
	// *********************************** <PERSISTENT VALUES> ***********************************

	NUM *g_lastM = g_last_lines + 3 * MAX_HAPLEN * blockIdx.x;
	NUM *g_lastX = g_lastM + MAX_HAPLEN;
	NUM *g_lastY = g_lastX + MAX_HAPLEN;

	__shared__ NUM s_xp[nThreads], s_yp[nThreads], s_mp[nThreads];
	__shared__ NUM s_xpp[nThreads], s_ypp[nThreads], s_mpp[nThreads];

	// *********************************** </PERSISTENT VALUES> **********************************
	// *******************************************************************************************

	for (;;)
	{
		// *******************************************************************************************
		// ************************************** <PICK A PAIR> **************************************
		__syncthreads();

		__shared__ Pair s_pair;
		__shared__ int s_pair_index;

		if (threadIdx.x == 0)
		{
			s_pair_index = atomicAdd(g_pair_index, 1);
			if (s_pair_index < num_of_pairs)
				s_pair = g_pair[s_pair_index];
		}

		__syncthreads();

		int pair_index = s_pair_index;

		if (pair_index >= num_of_pairs)
			break;

		Pair pair = s_pair;

		if (pair.status == DONE)
			continue;

		// ************************************** </PICK A PAIR> *************************************
		// *******************************************************************************************

		int n_groups_of_rows = ((pair.rslen + 1) + (nThreads - 1)) / nThreads;
		for (int group_of_rows = 0; group_of_rows < n_groups_of_rows; group_of_rows++)
		{
			__syncthreads();

			int row = group_of_rows * nThreads + threadIdx.x;

			/******************************************************************************************
			*********** <SET THE VALUES THAT ARE CONSTANT DURING THE CALCULATION OF THE ROW> *********/
			char rs;
			char4 qidc;
			NUM mm, gm, mx, xx, my, yy, pq;
			if (row > 0 && row <= pair.rslen)
			{
				rs = (g_chunk + pair.offset_rs)[row-1]; // GLOBAL MEMORY ACCESS
				qidc = reinterpret_cast<char4 *>(g_chunk + pair.offset_qidc)[row-1]; // GLOBAL MEMORY ACCESS //?? each byte in the qidc char4 should be masked with 127 during chunk creation
				mm = NUM(1.0) - PH2PR(qidc.y) * PH2PR(qidc.z);
				gm = NUM(1.0) - PH2PR(qidc.w);
				mx = PH2PR(qidc.y);
				xx = PH2PR(qidc.w);
				my = (row == pair.rslen) ? NUM(1.0) : PH2PR(qidc.z);
				yy = (row == pair.rslen) ? NUM(1.0) : PH2PR(qidc.w);
				pq = PH2PR(qidc.x);
			}

			NUM k = INITIAL_CONSTANT<NUM>() / pair.haplen;

			/********** </SET THE VALUES THAT ARE CONSTANT DURING THE CALCULATION OF THE ROW> **********
			*******************************************************************************************/

			// diagonal 0
			s_mpp[threadIdx.x] = NUM(0.0);
			s_xpp[threadIdx.x] = NUM(0.0);
			s_ypp[threadIdx.x] = NUM(0.0);

			//diagonal 1
			char hap;
			NUM m, x, y, sum_m_x = NUM(0.0), coef;

			s_mp[threadIdx.x] = NUM(0.0);
			s_xp[threadIdx.x] = NUM(0.0);
			s_yp[threadIdx.x] = NUM(0.0);

			__syncthreads();

			/*
				The purpose of the code inside the following "if" is:

				* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = *
				* to properly set the values of the first two diagonals. There two situations *
				* that we need to take care:                                                  *
				* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = *

				a) the group of rows (i.e.: the horizontal _tile_) is the first one
				b) the group of rows is not the first one.

				In case (a), the values of x, y and m should be initialized with 0 for all
				rows but the row 0, in which m and x should be 0 and y should be equal to
				the constant "k" (defined inside the if).

				In case (b), the only cell with values different to 0 is the cell located
				at the first row of the group of rows (i.e.: the row traversed by the thread 0)
				and the column 1 (i.e.: the first cell of the diagonal 1).

				The values of x, y and m at that cell will depend on the last values computed
				in the previous group of rows.
			*/	

			if (threadIdx.x == 0) 
			{
				if (group_of_rows == 0) // case (a): 
				{
					s_ypp[0] = s_yp[0] = k;
				}
				else // case (b): just need to compute the values of m, x and y at [0][1]. NEED: the character [0] of the haplotype.
				{
					hap = g_chunk[pair.offset_hap]; // GLOBAL MEMORY ACCESS
					coef = pq;
					if (rs == hap || rs == 'N' || hap == 'N')
						coef = NUM(1.0) - coef;

					NUM m_up_left = g_lastM[0]; // lastM, lastX and lastY have (pair.haplen+1) cells, corresponding to all the columns of the matrix. The first column represents the empty string.
					NUM x_up_left = g_lastX[0];
					NUM y_up_left = g_lastY[0];
					NUM m_up = g_lastM[1]; 
					NUM x_up = g_lastX[1];

					m = coef * (m_up_left * mm + x_up_left * gm + y_up_left * gm);
					x = m_up * mx + x_up * xx;
					y = NUM(0.0);

					s_xp[0] = x;
					s_yp[0] = y;
					s_mp[0] = m;

					sum_m_x += m + x;
				}
			} // end of "if (threadIdx.x==0) 
		
			__syncthreads();

			int num_of_diagonals = (pair.haplen+1) + (pair.rslen+1) - 1;
			for (int d = 2; d < num_of_diagonals; d++) 
			{
				__syncthreads();

				int col = d - threadIdx.x;
				hap = 0;
				NUM m_up_left = NUM(0.0), x_up_left = NUM(0.0), y_up_left = NUM(0.0);
				NUM m_up = NUM(0.0), x_up = NUM(0.0);
				NUM m_left = s_mp[0], y_left = s_yp[0];

//__syncthreads();
				if (threadIdx.x==0)
				{
					if (group_of_rows == 0)
					{
						yy = NUM(1.0);
						my = NUM(0.0);
					}
					else
					{
						if (col > 0 && col <= pair.haplen)
						{
							hap = g_chunk[pair.offset_hap + (col-1)]; // GLOBAL MEMORY ACCESS
							m_up_left = g_lastM[col-1]; // GLOBAL MEMORY ACCESS
							x_up_left = g_lastX[col-1]; // GLOBAL MEMORY ACCESS
							y_up_left = g_lastY[col-1]; // GLOBAL MEMORY ACCESS
							m_up = g_lastM[col]; // GLOBAL MEMORY ACCESS
							x_up = g_lastX[col]; // GLOBAL MEMORY ACCESS
							m_left = s_mp[0];
							y_left = s_yp[0];
						}
					}
				} 
				else
				{
					if (col > 0 && col <= pair.haplen) // ?? haplotype should be padded and this "if", removed.
						hap = g_chunk[pair.offset_hap + (col-1)]; // GLOBAL MEMORY ACCESS
					m_up_left = s_mpp[threadIdx.x-1]; 
					x_up_left = s_xpp[threadIdx.x-1];
					y_up_left = s_ypp[threadIdx.x-1];
					m_up = s_mp[threadIdx.x-1];
					x_up = s_xp[threadIdx.x-1];
					m_left = s_mp[threadIdx.x];
					y_left = s_yp[threadIdx.x];
				}

				coef = pq;
				if (rs == hap || rs == 'N' || hap == 'N')
					coef = NUM(1.0) - coef;
				m = coef * (m_up_left * mm + x_up_left * gm + y_up_left * gm);
				x = m_up * mx + x_up * xx;
				y = m_left * my + y_left * yy;

				if ((threadIdx.x == (nThreads-1)) && (group_of_rows < n_groups_of_rows - 1) && (col >= 0) && (col <= pair.haplen))
				{
					g_lastX[col] = x; // GLOBAL MEMORY ACCESS
					g_lastY[col] = y; // GLOBAL MEMORY ACCESS
					g_lastM[col] = m; // GLOBAL MEMORY ACCESS
				}

				__syncthreads(); // ?? this __syncthreads() has shown to be necessary, but it is not clear to me the reason. TO DO: Understand deeply what's happening in here. This code works, but it is mandatory to know why removing this __syncthreads() makes the program fail.

				if (col >= 0 && col <= pair.haplen)
					sum_m_x += m + x;

				s_xpp[threadIdx.x] = s_xp[threadIdx.x];	
				s_ypp[threadIdx.x] = s_yp[threadIdx.x];	
				s_mpp[threadIdx.x] = s_mp[threadIdx.x];
				s_xp[threadIdx.x] = x; 
				s_yp[threadIdx.x] = y; 
				s_mp[threadIdx.x] = m;
			} // end of the for (d = 2 to {number of diagonals}). Always followed by __syncthreads();

			__syncthreads();

			if (row == pair.rslen) 
			{
				pair.result = double(log10(sum_m_x) - log10(INITIAL_CONSTANT<NUM>()));
				pair.status = (sum_m_x >= MIN_ACCEPTED<NUM>()) ? DONE : NOT_DONE;
				g_pair[s_pair_index] = pair; // GLOBAL MEMORY ACCESS
			}
		} // end of the for (group_of rows = ...). Always followed by __syncthreads();

		__syncthreads();

	} // end of the for (;;)

	return;
}

int create_chunk(vector<Pair> &pairs, string &chunk)
{
	vector<Pair>().swap(pairs); 
	chunk.clear();

	int current_offset = 0;
	std::string hap, rs, q, i, d, c, i1, i2, tchunk;
	while (pairs.size() < NUM_OF_TESTCASES_PER_ITERATION && (cin >> hap >> rs >> q >> i >> d >> c >> i1 >> i2).good())
	{
		tchunk.clear();
		tchunk = hap + string(1, '\0') + rs + string(1, '\0');
		tchunk += string(GAP_TO(tchunk.size()/*tchunk_sz*/, 4), '\0');

		for (int x = 0; x < rs.size(); x++)
		{
			char tq = q[x] - 33; if (tq < 6) tq = 6; tq = (tq & 127);
			char ti = (i[x] - 33) & 127;
			char td = (d[x] - 33) & 127;
			char tc = (c[x] - 33) & 127;
			tchunk += string(1, tq) + string(1, ti) + string(1, td) + string(1, tc);
		}
		tchunk += string(4, '\0');
		tchunk += string(GAP_TO(tchunk.size(), 128), '\0');

		Pair p;
		p.status = NOT_DONE;
		p.result = 0.0;
		p.offset_hap = current_offset;
		p.offset_rs = current_offset + (hap.size()+1);
		int sz_hap_rs = (hap.size()+1) + (rs.size()+1);
		p.offset_qidc = current_offset +  sz_hap_rs + GAP_TO(sz_hap_rs, 4);
		p.haplen = hap.size();
		p.rslen = rs.size();

		chunk += tchunk;
		pairs.push_back(p);

		current_offset += tchunk.size();
		hap.clear();
		rs.clear();
		q.clear();
		i.clear();
		d.clear();
		c.clear();
	}

	return pairs.size();
}

class Timing
{
public:
	Timing(string t) : st(now()), tot(0.0), title(t)
	{
	}

	~Timing()
	{
		std::cerr << title << tot << " seconds\n";
	}

	void start()
	{
		st = now();
	}

	void acc()
	{
		tot += (now() - st);
	}

private:
	static double now()
	{
		struct timeval v;
		gettimeofday(&v, (struct timezone *) NULL);
		return v.tv_sec + v.tv_usec/1.0e6;
	}

	double st, tot;
	string title;
};

int main(void)
{
	Timing TotalTime(string("TOTAL: "));
	Timing ComputationTime(string("COMPUTATION: "));

	dim3 gridDim(NBLOCKS);
	dim3 blockDim(NTHREADS);

	string chunk;
	vector<Pair> pairs;

	while (create_chunk(pairs, chunk))
	{
		ComputationTime.start();

		char *g_chunk;
		int padd_sz = 128;
		gpuErrchk( cudaMalloc(&g_chunk, chunk.size() + 2 * padd_sz) );

		gpuErrchk( cudaMemset(g_chunk, 0, padd_sz) );
		gpuErrchk( cudaMemset(g_chunk + padd_sz + chunk.size(), 0, padd_sz) );
		gpuErrchk( cudaMemcpy(g_chunk + padd_sz, chunk.c_str(), chunk.size(), cudaMemcpyHostToDevice) );
		g_chunk += padd_sz;

		Pair *g_pair;
		gpuErrchk( cudaMalloc(&g_pair, pairs.size() * sizeof(Pair)) );
		gpuErrchk( cudaMemcpy(g_pair, &(pairs[0]), pairs.size() * sizeof(Pair), cudaMemcpyHostToDevice) );

		int *g_pair_index;
		gpuErrchk( cudaMalloc(&g_pair_index, sizeof(int)) );

//cudaPrintfInit();

		void *g_last_lines;
		gpuErrchk( cudaMalloc(&g_last_lines, NBLOCKS * MAX_HAPLEN * 3 * sizeof(double)) );

		gpuErrchk( cudaMemset(g_pair_index, 0, sizeof(int)) );
		compute_full_scores<NTHREADS, float><<<gridDim, blockDim>>>(g_chunk, pairs.size(), g_pair, g_pair_index, reinterpret_cast<float *>(g_last_lines));

cudaThreadSynchronize(); // agregado por si acaso... no es esto...
		gpuErrchk( cudaMemset(g_pair_index, 0, sizeof(int)) );
		compute_full_scores<NTHREADS, double><<<gridDim, blockDim>>>(g_chunk, pairs.size(), g_pair, g_pair_index, reinterpret_cast<double *>(g_last_lines));

//cudaPrintfDisplay();
//cudaPrintfEnd();

		gpuErrchk( cudaMemcpy(&(pairs[0]), g_pair, pairs.size() * sizeof(Pair), cudaMemcpyDeviceToHost) );

		cout << std::setprecision(16);
		for (int p = 0; p < pairs.size(); p++)
			cout << pairs[p].result << "\n";

		gpuErrchk( cudaFree(g_last_lines) );
		gpuErrchk( cudaFree(g_pair) );
		gpuErrchk( cudaFree(g_pair_index) );
		gpuErrchk( cudaFree(g_chunk-padd_sz) );

		ComputationTime.acc();
	}

	TotalTime.acc();

	return 0;
}

