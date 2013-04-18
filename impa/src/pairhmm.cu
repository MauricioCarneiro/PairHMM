#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include "common.h"

typedef unsigned long ul;
typedef unsigned char uc;
typedef double dbl;

#define BLOCKWIDTH 64
#define BLOCKHEIGHT 1
#define COMPDIAGS 60
#define COMPBUFFSIZE 30

#define MIN(a,b) (((a)<(b))?(a):(b))
#define PARALLELMXY (BLOCKHEIGHT == 3)
#define TX (threadIdx.x)
#define TY (threadIdx.y)
#define TID (threadIdx.x + BLOCKWIDTH * threadIdx.y)
#define NTH (BLOCKHEIGHT * BLOCKWIDTH)
#define THREADM (threadIdx.y == 0)
#define THREADX (threadIdx.y == 1)
#define THREADY (threadIdx.y == 2)

#define MAX_COMPARISONS_PER_SPLIT ((unsigned long) 35000)
#define MAX_H 5120
#define MAX_R 1536
#define MAX_SIMULTANEOUS_BLOCKS 150

template<ul ROWS, ul DIAGS, ul BUFFSIZE>
struct PubVars
{
	ul H, R;
	char h[ROWS+DIAGS-1];
	char r[ROWS];
	uc q[ROWS], i[ROWS], d[ROWS], c[ROWS];
	ul nblockrows, nblockcols;
	Haplotype hap;
	ReadSequence rs;
	char *chunk;
	dbl ph2pr[128];
	dbl m1[ROWS], m2[ROWS], m3[ROWS];
	dbl x1[ROWS], x2[ROWS], x3[ROWS];
	dbl y1[ROWS], y2[ROWS], y3[ROWS];
	dbl *m, *mp, *mpp;
	dbl *x, *xp, *xpp;
	dbl *y, *yp, *ypp;
	dbl *g_lastM, *g_lastX, *g_lastY;
	dbl lastM[DIAGS+1], lastX[DIAGS+1], lastY[DIAGS+1];
	dbl buffM[BUFFSIZE], buffX[BUFFSIZE], buffY[BUFFSIZE];
	ul buffsz;
	ul buffstart;
};

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void flush_buffer(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	for (int k = TID; k < pv.buffsz; k += NTH)
	{
		pv.g_lastM[pv.buffstart + k] = pv.buffM[k];
		pv.g_lastX[pv.buffstart + k] = pv.buffX[k];
		pv.g_lastY[pv.buffstart + k] = pv.buffY[k];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void load_previous_results(int &j, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	for (int c = ((j == 0) ? 1 : 0) + TID; c < MIN(pv.H + 1 - j * DIAGS, DIAGS) + 1; c += NTH)
	{
		pv.lastM[c] = pv.g_lastM[j * DIAGS + c - 1];
		pv.lastX[c] = pv.g_lastX[j * DIAGS + c - 1];
		pv.lastY[c] = pv.g_lastY[j * DIAGS + c - 1];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void load_hap_data(int &j, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	if (TID < DIAGS)
		for (int c = TID; c < ROWS - 1; c += DIAGS)
			pv.h[c] = pv.h[c + DIAGS];

	__syncthreads();

	for (int c = ((j == 0) ? 1 : 0) + TID; c < MIN(DIAGS, pv.H + 1 - j * DIAGS); c += NTH)
		pv.h[ROWS - 1 + c] = pv.chunk[pv.hap.h + j * DIAGS + c - 1];
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void load_read_data(int &i, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	for (int r = ((i == 0) ? 1 : 0) + TID; r < MIN(ROWS, pv.R + 1 - i * ROWS); r += NTH)
	{
		pv.r[r] = pv.chunk[pv.rs.r + i * ROWS + r - 1];
		pv.q[r] = pv.chunk[pv.rs.qual + i * ROWS + r - 1];
		pv.i[r] = pv.chunk[pv.rs.ins + i * ROWS + r - 1];
		pv.d[r] = pv.chunk[pv.rs.del + i * ROWS + r - 1];
		pv.c[r] = pv.chunk[pv.rs.cont + i * ROWS + r - 1];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void notfirstline_firstcolum(int &r, PubVars<ROWS, DIAGS, BUFFSIZE> &pv, dbl &_m, dbl &_x)
{
	if (PARALLELMXY)
	{
		if (THREADM)
			pv.m[r] = 0.0;

		if (THREADX)
			pv.x[r] = _m * pv.ph2pr[pv.i[r]] + _x * pv.ph2pr[pv.c[r]];

		if (THREADY)
			pv.y[r] = 0.0;
	}
	else
	{
		pv.m[r] = 0.0;
		pv.x[r] = _m * pv.ph2pr[pv.i[r]] + _x * pv.ph2pr[pv.c[r]];
		pv.y[r] = 0.0;
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void notfirstline_notfirstcolumn(int &r, int &i, int &diag, PubVars<ROWS, DIAGS, BUFFSIZE> &pv, dbl &M_m, dbl &M_x, dbl &M_y, dbl &X_m, dbl &X_x, dbl &Y_m, dbl &Y_y)
{
	dbl t1, t2, t3, dist;

	if (PARALLELMXY)
	{
		if (THREADM)
		{
			t1 = 1.0 - pv.ph2pr[(pv.i[r] + pv.d[r]) & 127];
			t2 = 1.0 - pv.ph2pr[pv.c[r]];
			t3 = pv.ph2pr[pv.q[r]];
			dist = (pv.r[r] == pv.h[ROWS - 1 + diag - r] || pv.r[r] == 'N' || pv.h[ROWS - 1 + diag - r] == 'N') ? 1.0 - t3 : t3;
			pv.m[r] = dist * (M_m * t1 + M_x * t2 + M_y * t2);
		}

		if (THREADX)
		{
			t1 = pv.ph2pr[pv.i[r]];
			t2 = pv.ph2pr[pv.c[r]];
			pv.x[r] = X_m * t1 + X_x * t2;
		}

		if (THREADY)
		{
			t1 = ((unsigned)r + i * ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.d[r]];
			t2 = ((unsigned)r + i * ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.c[r]];
			pv.y[r] = Y_m * t1 + Y_y * t2;
		}
	}
	else
	{
		t1 = 1.0 - pv.ph2pr[(pv.i[r] + pv.d[r]) & 127];
		t2 = 1.0 - pv.ph2pr[pv.c[r]];
		t3 = pv.ph2pr[pv.q[r]];
		dist = (pv.r[r] == pv.h[ROWS - 1 + diag - r] || pv.r[r] == 'N' || pv.h[ROWS - 1 + diag - r] == 'N') ? 1.0 - t3 : t3;
		pv.m[r] = dist * (M_m * t1 + M_x * t2 + M_y * t2);

		t1 = pv.ph2pr[pv.i[r]];
		t2 = pv.ph2pr[pv.c[r]];
		pv.x[r] = X_m * t1 + X_x * t2;

		t1 = ((unsigned)r + i * ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.d[r]];
		t2 = ((unsigned)r + i * ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.c[r]];
		pv.y[r] = Y_m * t1 + Y_y * t2;
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void firstline_firstcolum(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	if (PARALLELMXY)
	{
		if (THREADM)
			pv.m[0] = 1.0e300 / (pv.H > pv.R ? pv.H - pv.R + 1 : 1);

		if (THREADX)
			pv.x[0] = 0.0;

		if (THREADY)
			pv.y[0] = 0.0;
	}
	else
	{
		pv.m[0] = 1.0e300 / (pv.H > pv.R ? pv.H - pv.R + 1 : 1);
		pv.x[0] = 0.0;
		pv.y[0] = 0.0;
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void firstline_notfirstcolum(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	if (PARALLELMXY)
	{
		if (THREADM)
			pv.m[0] = 0.0;

		if (THREADX)
			pv.x[0] = 0.0;

		if (THREADY)
			pv.y[0] = pv.mp[0] + pv.yp[0];
	}
	else
	{
		pv.m[0] = 0.0;
		pv.x[0] = 0.0;
		pv.y[0] = pv.mp[0] + pv.yp[0];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void rotatediags(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	dbl *sw;
	sw = pv.mpp; pv.mpp = pv.mp; pv.mp = pv.m; pv.m = sw;
	sw = pv.xpp; pv.xpp = pv.xp; pv.xp = pv.x; pv.x = sw;
	sw = pv.ypp; pv.ypp = pv.yp; pv.yp = pv.y; pv.y = sw;
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__device__ inline void block(int &i, int &j, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	if (i > 0)
		load_previous_results<ROWS, DIAGS, BUFFSIZE>(j, pv);

	load_hap_data<ROWS, DIAGS, BUFFSIZE>(j, pv);

	__syncthreads();

	int nRows = MIN(ROWS, (int)(pv.R + 1 - i * ROWS));
	for (int diag = 0; diag < DIAGS; diag++)
	{
		int r, c;
		for (r = TX; r < nRows; r += BLOCKWIDTH)
		{
			c = diag - r + j * DIAGS;
			if (c >= 0 && c <= (int)pv.H)
				if (r == 0)
					if (i == 0)
						if (c == 0)
							firstline_firstcolum<ROWS, DIAGS, BUFFSIZE>(pv);
						else
							firstline_notfirstcolum<ROWS, DIAGS, BUFFSIZE>(pv);
					else
						if (c == 0)
							notfirstline_firstcolum<ROWS, DIAGS, BUFFSIZE>(r, pv, pv.lastM[1], pv.lastX[1]);
						else
							notfirstline_notfirstcolumn<ROWS, DIAGS, BUFFSIZE>(r, i, diag, pv, pv.lastM[diag], pv.lastX[diag], pv.lastY[diag], pv.lastM[diag + 1], pv.lastX[diag + 1], pv.mp[r], pv.yp[r]);
				else
					if (c == 0)
						notfirstline_firstcolum<ROWS, DIAGS, BUFFSIZE>(r, pv, pv.mp[r-1], pv.xp[r-1]);
					else
						notfirstline_notfirstcolumn<ROWS, DIAGS, BUFFSIZE>(r, i, diag, pv, pv.mpp[r-1], pv.xpp[r-1], pv.ypp[r-1], pv.mp[r-1], pv.xp[r-1], pv.mp[r], pv.yp[r]);
		}

		__syncthreads();

		r = nRows - 1;
		c = diag - r + j * DIAGS;
		if ((TID == 0) && (c >= 0) && (c <= (int)pv.H) && (i < pv.nblockrows - 1))
		{
			pv.buffM[pv.buffsz] = pv.m[r];
			pv.buffX[pv.buffsz] = pv.x[r];
			pv.buffY[pv.buffsz] = pv.y[r];
			pv.buffsz++;
		}

		__syncthreads();

		if (pv.buffsz == BUFFSIZE)
		{
			flush_buffer<ROWS, DIAGS, BUFFSIZE>(pv);

			if (TID == 0)
			{
				pv.buffstart += BUFFSIZE;
				pv.buffsz = 0;
			}
			__syncthreads();
		}

		if (TID == 0)
			rotatediags<ROWS, DIAGS, BUFFSIZE>(pv);

		__syncthreads();
	}

	return;
}

template <int ROWS, int DIAGS, int BUFFSIZE>
__global__ void compare(Memory mem, dbl *g_lastLinesArr, int *g_lastLinesIndex, int *g_compIndex)
{
	__shared__ PubVars<ROWS, DIAGS, BUFFSIZE> pv;
	__shared__ int compIndex;

	for (int i = TID; i < 128; i += NTH)
		pv.ph2pr[i] = pow(10.0, -((double) i) / 10.0);

	if (TID == 0)
	{
		int lastLinesIndex = atomicAdd(g_lastLinesIndex, 1);
		pv.g_lastM = g_lastLinesArr + (lastLinesIndex * 3 * MAX_H);
		pv.g_lastX = pv.g_lastM + MAX_H;
		pv.g_lastY = pv.g_lastX + MAX_H;
		pv.chunk = mem.chunk;
	}

	__syncthreads();

	for (;;)
	{
		if (TID == 0)
			compIndex = atomicAdd(g_compIndex, 1);

		__syncthreads();

		if (compIndex >= mem.nres)
			break;

		if (TID == 0)
		{
			pv.rs = mem.r[mem.cmpR[compIndex]];
			pv.hap = mem.h[mem.cmpH[compIndex]];
			pv.buffsz = 0;
			pv.buffstart = 0;
			pv.H = pv.hap.H;
			pv.R = pv.rs.R;
			pv.nblockrows = (pv.R + ROWS) / ROWS;
			pv.nblockcols = (ROWS + pv.H + DIAGS - 1) / DIAGS;
			pv.m = pv.m1; pv.mp = pv.m2; pv.mpp = pv.m3;
			pv.y = pv.y1; pv.yp = pv.y2; pv.ypp = pv.y3;
			pv.x = pv.x1; pv.xp = pv.x2; pv.xpp = pv.x3;
		}

		__syncthreads();

		for (int i = 0; i < pv.nblockrows; i++)
		{
			load_read_data<ROWS, DIAGS, BUFFSIZE>(i, pv);

			if (TID == 0)
			{
				pv.buffstart = 0;
				pv.buffsz = 0;
			}

			__syncthreads();

			for (int j = 0; j < pv.nblockcols; j++)
			{
				block<ROWS, DIAGS, BUFFSIZE>(i, j, pv);
				__syncthreads();
			}
			flush_buffer<ROWS, DIAGS, BUFFSIZE>(pv);

			__syncthreads();
		}

		if (TID == 0)
		{
			int resr = (int)pv.R - (int)(pv.nblockrows - 1) * ROWS;
			int swaps = (DIAGS * (int)pv.nblockcols - resr - (int)(pv.H + 1)) % 3;

			if (swaps == 2)
				mem.res[compIndex] = log10(pv.m[resr] + pv.x[resr] + pv.y[resr]) - 300.0;
			else if (swaps == 0)
				mem.res[compIndex] = log10(pv.mp[resr] + pv.xp[resr] + pv.yp[resr]) - 300.0;
			else
				mem.res[compIndex] = log10(pv.mpp[resr] + pv.xpp[resr] + pv.ypp[resr]) - 300.0;
		}

		__syncthreads();
	}

	return;
}

int split(Memory &h_big, Memory &ret)
{
	static ul lastGroup = 0;
	static ul offset_res = 0;
	ul offset_h = h_big.g[lastGroup].fstH;
	ul offset_r = h_big.g[lastGroup].fstR;
	ul chunk_begin, chunk_end;
	ul j;
	ul fstG = lastGroup;

	if (lastGroup >= h_big.ng)
		return 1;

	ret.nres = 0;
	ret.ng = 0;
	ret.nh = 0;
	ret.nr = 0;
	ret.chunk_sz = 0;

	while ((lastGroup < h_big.ng) && (ret.nres + h_big.g[lastGroup].nR * h_big.g[lastGroup].nH < MAX_COMPARISONS_PER_SPLIT))
	{
		ret.nres += h_big.g[lastGroup].nR * h_big.g[lastGroup].nH;
		ret.ng++;
		ret.nh += h_big.g[lastGroup].nH;
		ret.nr += h_big.g[lastGroup].nR;
		lastGroup++;
	}

	if (ret.nres == 0)
	{
		fprintf(stderr, "There exists a group with more than MAX_COMPARISONS_PER_SPLIT comparisons\n");
		exit(0);
	}

	chunk_begin = h_big.r[h_big.g[fstG].fstR].r;
	chunk_end = h_big.h[h_big.g[lastGroup-1].fstH + h_big.g[lastGroup-1].nH-1].h + h_big.h[h_big.g[lastGroup-1].fstH + h_big.g[lastGroup-1].nH-1].H + 1;
	ret.chunk_sz = (chunk_end - chunk_begin + 1);
	ret.chunk = h_big.chunk + chunk_begin;

	for (j = 0; j < ret.nh; j++)
	{
		ret.h[j] = h_big.h[offset_h + j];
		ret.h[j].h -= chunk_begin;
	}

	for (j = 0; j < ret.nr; j++)
	{
		ret.r[j] = h_big.r[offset_r + j];
		ret.r[j].r -= chunk_begin;
		ret.r[j].qual -= chunk_begin;
		ret.r[j].ins -= chunk_begin;
		ret.r[j].del -= chunk_begin;
		ret.r[j].cont -= chunk_begin;
	}

	for (j = 0; j < ret.nres; j++)
	{
		ret.cmpH[j] = h_big.cmpH[offset_res + j] - offset_h;
		ret.cmpR[j] = h_big.cmpR[offset_res + j] - offset_r;
	}

	offset_res += ret.nres;

	return 0;
}

int main(int argc, char **argv)
{
	Memory h_big, h_small, d_mem;
	ul already = 0;
	dbl *g_lastlines;
	int *g_compIndex, compIndex = 0;
	int *g_lastLinesIndex, lastLinesIndex = 0;

	struct
	{
		double start, init_memory, mallocs, comp, output, end;
		double t1, t2;
		float kernel;
	} times;
	times.kernel = 0.f;

	times.start = right_now();
	times.t1 = right_now();
	init_memory(argv[1], &h_big);
	times.t2 = right_now();
	times.init_memory = times.t2 - times.t1;

	times.t1 = right_now();

	h_small.r = (ReadSequence *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(ReadSequence));
	h_small.h = (Haplotype *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(Haplotype));
	h_small.chunk = (char *) malloc(MAX_COMPARISONS_PER_SPLIT * (MAX_H + MAX_R * 5));
	h_small.cmpH = (ul *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	h_small.cmpR = (ul *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	h_small.res = (dbl *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(dbl));
	h_small.g = NULL;
	h_small.phred_to_prob = NULL;

	g_lastlines = NULL;
	g_compIndex = NULL;
	g_lastLinesIndex = NULL;
	d_mem.r = NULL;
	d_mem.h = NULL;
	d_mem.chunk = NULL;
	d_mem.cmpH = NULL;
	d_mem.cmpR = NULL;
	d_mem.res = NULL;

	cudaMalloc(&g_lastlines, 3 * sizeof(dbl) * MAX_H * MAX_SIMULTANEOUS_BLOCKS);
	cudaMalloc(&g_compIndex, sizeof(int));
	cudaMalloc(&g_lastLinesIndex, sizeof(int));
	cudaMalloc(&(d_mem.r), MAX_COMPARISONS_PER_SPLIT * sizeof(ReadSequence));
	cudaMalloc(&(d_mem.h), MAX_COMPARISONS_PER_SPLIT * sizeof(Haplotype));
	cudaMalloc(&(d_mem.chunk), MAX_COMPARISONS_PER_SPLIT * (MAX_H + MAX_R * 5));
	cudaMalloc(&(d_mem.cmpH), MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	cudaMalloc(&(d_mem.cmpR), MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	cudaMalloc(&(d_mem.res), MAX_COMPARISONS_PER_SPLIT * sizeof(dbl));
	d_mem.g = NULL;
	d_mem.phred_to_prob = NULL;

	if (!g_lastLinesIndex || !g_lastlines || !g_compIndex || !d_mem.r || !d_mem.h || !d_mem.chunk || !d_mem.cmpH || !d_mem.cmpR || !d_mem.res)
	{
		fprintf(stderr, "Some malloc went wrong...\n");
		exit(0);
	}

	times.t2 = right_now();
	times.mallocs = times.t2 - times.t1;

	times.t1 = right_now();
	while (!split(h_big, h_small))
	{
		cudaEvent_t kernel_start, kernel_stop;
		float k_time;
		cudaEventCreate(&kernel_start);
		cudaEventCreate(&kernel_stop);

		d_mem.ng = h_small.ng;
		d_mem.nh = h_small.nh;
		d_mem.nr = h_small.nr;
		d_mem.chunk_sz = h_small.chunk_sz;
		d_mem.nres = h_small.nres;
		cudaMemcpy(g_compIndex, &compIndex, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(g_lastLinesIndex, &lastLinesIndex, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.r, h_small.r, h_small.nr * sizeof(ReadSequence), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.h, h_small.h, h_small.nh * sizeof(Haplotype), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.chunk, h_small.chunk, h_small.chunk_sz, cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.cmpH, h_small.cmpH, h_small.nres * sizeof(ul), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.cmpR, h_small.cmpR, h_small.nres * sizeof(ul), cudaMemcpyHostToDevice);
		dim3 gridDim(MAX_SIMULTANEOUS_BLOCKS);
		dim3 blockDim(BLOCKWIDTH, BLOCKHEIGHT);
		cudaEventRecord(kernel_start, 0);
		compare<BLOCKWIDTH, COMPDIAGS, COMPBUFFSIZE><<<gridDim, blockDim>>>(d_mem, g_lastlines, g_lastLinesIndex, g_compIndex);
		cudaEventRecord(kernel_stop, 0);
		cudaEventSynchronize(kernel_stop);

		cudaMemcpy(h_big.res + already, d_mem.res, d_mem.nres * sizeof(dbl), cudaMemcpyDeviceToHost);
		already += d_mem.nres;
		cudaEventElapsedTime(&k_time, kernel_start, kernel_stop);
		times.kernel += k_time;
	}
	times.t2 = right_now();
	times.comp = times.t2 - times.t1;

	cudaFree(g_lastlines);
	cudaFree(g_lastLinesIndex);
	cudaFree(g_compIndex);
	cudaFree(d_mem.r);
	cudaFree(d_mem.h);
	cudaFree(d_mem.chunk);
	cudaFree(d_mem.cmpH);
	cudaFree(d_mem.cmpR);
	cudaFree(d_mem.res);

	times.t1 = right_now();
	output(h_big.res, h_big.nres, argv[2]);
	times.t2 = right_now();
	times.output = times.t2 - times.t1;
	times.end = right_now();

	printf("INIT_MEMORY: %g\n", times.init_memory * 1000.0);
	printf("MALLOCS: %g\n", times.mallocs * 1000.0);
	printf("COMPUTATION: %g\n", times.comp * 1000.0);
	printf("KERNEL: %f\n", times.kernel);
	printf("OUTPUT: %g\n", times.output * 1000.0);
	printf("TOTAL (measured inside program): %g\n", (times.end - times.start) * 1000.0);

	return 0;
}

