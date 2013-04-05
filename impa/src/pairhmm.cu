#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <cctype>

using std::cout;
using std::isprint;

#include "common.h"
//#include "cuPrintf.cu"

typedef unsigned long ul;
typedef unsigned char uc;
typedef double dbl;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX_COMPARISONS_PER_SPLIT ((unsigned long) 250000)

#define MAX_H 10240
#define BLOCKS_PER_SM 8
#define SMS 14

#define TID (threadIdx.x + blockDim.x * threadIdx.y)
#define NTH (blockDim.x * blockDim.y)


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

struct Timers
{
	Timers()
	{
		split = 0;
		compute = 0;
		HtoD = 0;
		DtoH = 0;
	}

	void acc()
	{
		split += (t2 - t1);
		HtoD += (t3 - t2);
		compute += (t4 - t3);
		DtoH += (t5 - t4);
	}

	dbl start, t0, t1, t2, t3, t4, t5, t6, t7, split, compute, HtoD, DtoH;
};

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
__device__ void load_previous_results(ul j, PubVars<ROWS, DIAGS, BUFFSIZE> *pv)
{
	ul tid = threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = blockDim.x * blockDim.y;
	ul c;

	__syncthreads();
	for (c = ((j==0)?1:0) + tid; c < MIN(pv->H+1-j*DIAGS, DIAGS)+1; c+=nth)
	{
		pv->lastM[c] = pv->g_lastM[j*DIAGS+c-1];
		pv->lastX[c] = pv->g_lastX[j*DIAGS+c-1];
		pv->lastY[c] = pv->g_lastY[j*DIAGS+c-1];
	}
	__syncthreads();
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
__device__ void load_hap_data(ul j, PubVars<ROWS, DIAGS, BUFFSIZE> *pv)
{
	ul tid = threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = blockDim.x * blockDim.y;
	ul c;

	__syncthreads();
	for (c = tid; c < ROWS-1; c+=nth)
		pv->h[c] = pv->h[c + DIAGS];
	__syncthreads();
	for (c=((j==0)?1:0)+tid; c < MIN(DIAGS, pv->H + 1 - j * DIAGS); c+=nth)
		pv->h[ROWS-1+c] = pv->chunk[pv->hap.h + j*DIAGS + c - 1];
	__syncthreads();
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
__device__ void load_read_data(ul i, PubVars<ROWS, DIAGS, BUFFSIZE> *pv)
{
	ul tid = threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = blockDim.x * blockDim.y;
	ul r;

	__syncthreads();
	for (r = ((i==0)?1:0) + tid; r < MIN(ROWS, pv->R+1-i*ROWS); r+=nth)
	{
		pv->r[r] = pv->chunk[pv->rs.r + i*ROWS + r - 1];
		pv->q[r] = pv->chunk[pv->rs.qual + i*ROWS + r - 1];
		pv->i[r] = pv->chunk[pv->rs.ins + i*ROWS + r - 1];
		pv->d[r] = pv->chunk[pv->rs.del + i*ROWS + r - 1];
		pv->c[r] = pv->chunk[pv->rs.cont + i*ROWS + r - 1];
	}
	__syncthreads();
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
__device__ void block(ul i, ul j, PubVars<ROWS, DIAGS, BUFFSIZE> *pv)
{
	ul tid = threadIdx.x + blockDim.x * threadIdx.y;
	ul tx = threadIdx.x;
	ul ty = threadIdx.y;
	ul nth = blockDim.x * blockDim.y;

	dbl *sw;
	int c, r, nRows;
	ul diag, k;
	dbl _m, _x, _y, dist, t1, t2, t3;

	if (i > 0)
	{
		load_previous_results<ROWS, DIAGS, BUFFSIZE>(j, pv);
	}
	__syncthreads();

	load_hap_data<ROWS, DIAGS, BUFFSIZE>(j, pv);
	__syncthreads();

	nRows = MIN((int)ROWS, (int)(pv->R+1-i*ROWS));

	for (diag=0; diag < DIAGS; diag++)
	{
		for (r=tx; r < nRows; r+=blockDim.x)
		{
			c = diag - r + j * DIAGS;
			if (c < 0 || c > (int)pv->H)
				continue;
			if (r == 0)
			{
				if (i == 0)
				{
					if (c == 0)
					{
						if (ty == 0)
						{
							pv->m[r] = 1.0e300 / (pv->H > pv->R ? pv->H - pv->R + 1 : 1);
						}
						else if (ty == 1)
						{
							pv->x[r] = 0.0;
						}
						else
						{
							pv->y[r] = 0.0;
						}
					}
					else
					{
						if (ty == 0)
						{
							pv->m[r] = 0.0;
						}
						else if (ty == 1)
						{
							pv->x[r] = 0.0;
						}
						else
						{
							pv->y[r] = pv->mp[0] + pv->yp[0];
						}
					}
				}
				else
				{
					if (c == 0)
					{
						if (ty == 0)
						{
							pv->m[r] = 0.0;
						}
						else if (ty == 1)
						{
							t1 = pv->ph2pr[pv->i[r]];
							t2 = pv->ph2pr[pv->c[r]];
							_m = pv->lastM[1];
							_x = pv->lastX[1];
							pv->x[r] = _m * t1 + _x * t2;
						}
						else
						{
							pv->y[r] = 0.0;
						}
					}
					else
					{
						if (ty == 0)
						{
							t1 = 1.0 - pv->ph2pr[(pv->i[r] + pv->d[r]) & 127];
							t2 = 1.0 - pv->ph2pr[pv->c[r]];
							t3 = pv->ph2pr[pv->q[r]];
							_m = pv->lastM[diag];
							_x = pv->lastX[diag];
							_y = pv->lastY[diag];
							dist = (pv->r[r]==pv->h[ROWS-1+diag-r] || pv->r[r]=='N' || pv->h[ROWS-1+diag-r]=='N') ? 1.0-t3 : t3;
							pv->m[r] = dist * (_m * t1 + _x * t2 + _y * t2);
						}
						else if (ty == 1)
						{
							t1 = pv->ph2pr[pv->i[r]];
							t2 = pv->ph2pr[pv->c[r]];
							_m = pv->lastM[diag+1];
							_x = pv->lastX[diag+1];
							pv->x[r] = _m * t1 + _x * t2;
						}
						else
						{
							t1 = ((unsigned)r+i*ROWS == pv->R) ? 1.0 : pv->ph2pr[pv->d[r]];
							t2 = ((unsigned)r+i*ROWS == pv->R) ? 1.0 : pv->ph2pr[pv->c[r]];
							_m = pv->mp[r];
							_y = pv->yp[r];
							pv->y[r] = _m * t1 + _y * t2;
						}
					}
				}
			}
			else
			{

				if (c == 0)
				{
					if (ty == 0)
					{
						pv->m[r] = 0.0;
					}
					else if (ty == 1)
					{
						t1 = pv->ph2pr[pv->i[r]];
						t2 = pv->ph2pr[pv->c[r]];
						_m = pv->mp[r-1];
						_x = pv->xp[r-1];
						pv->x[r] = _m * t1 + _x * t2;
					}
					else
					{
						pv->y[r] = 0.0;
					}
				}
				else
				{
					if (ty == 0)
					{
						t1 = 1.0 - pv->ph2pr[(pv->i[r] + pv->d[r]) & 127];
						t2 = 1.0 - pv->ph2pr[pv->c[r]];
						t3 = pv->ph2pr[pv->q[r]];
						_m = pv->mpp[r-1];
						_x = pv->xpp[r-1];
						_y = pv->ypp[r-1];
						dist = (pv->r[r]==pv->h[ROWS-1+diag-r] || pv->r[r]=='N' || pv->h[ROWS-1+diag-r]=='N') ? 1.0-t3 : t3;
						pv->m[r] = dist * (_m * t1 + _x * t2 + _y * t2);
					}
					else if (ty == 1)
					{
						t1 = pv->ph2pr[pv->i[r]];
						t2 = pv->ph2pr[pv->c[r]];
						_m = pv->mp[r-1];
						_x = pv->xp[r-1];
						pv->x[r] = _m * t1 + _x * t2;
					}
					else
					{
						t1 = ((unsigned)r+i*ROWS == pv->R) ? 1.0 : pv->ph2pr[pv->d[r]];
						t2 = ((unsigned)r+i*ROWS == pv->R) ? 1.0 : pv->ph2pr[pv->c[r]];
						_m = pv->mp[r];
						_y = pv->yp[r];
						pv->y[r] = _m * t1 + _y * t2;
					}
				}
			}
			if ((r == nRows-1) && (i < pv->nblockrows-1) && (ty == 0))
			{
				(pv->buffM)[pv->buffsz] = pv->m[r];
				(pv->buffX)[pv->buffsz] = pv->x[r];
				(pv->buffY)[pv->buffsz] = pv->y[r];
				(pv->buffsz)++;
			}
			__syncthreads();
		}
		__syncthreads();
		if (pv->buffsz==BUFFSIZE)
		{
			for (k=tid; k < pv->buffsz; k+=nth)
			{
				pv->g_lastM[pv->buffstart + k] = pv->buffM[k];
				pv->g_lastX[pv->buffstart + k] = pv->buffX[k];
				pv->g_lastY[pv->buffstart + k] = pv->buffY[k];
			}
			if (tid == 0)
			{
				pv->buffstart+=BUFFSIZE;
				pv->buffsz=0;
			}
		}

		__syncthreads();
		if (tid == 0)
		{
			sw=pv->mpp; pv->mpp=pv->mp; pv->mp=pv->m; pv->m=sw;
			sw=pv->xpp; pv->xpp=pv->xp; pv->xp=pv->x; pv->x=sw;
			sw=pv->ypp; pv->ypp=pv->yp; pv->yp=pv->y; pv->y=sw;
		}
		__syncthreads();
	}

	return;
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
__global__ void compare(Memory mem, dbl *g_lines, int *g_counter)
{
	ul i, j, k;
	__shared__ PubVars<ROWS, DIAGS, BUFFSIZE> pv;
	__shared__ ul compIndex;

	////////////////////////   ONLY ONCE   ///////////////////////////////
	for (i=0+TID; i < 128; i+=NTH)
		pv.ph2pr[i] = pow(10.0, -((double) i) / 10.0);

	if (TID == 0)
	{
		pv.g_lastM = g_lines + (blockIdx.x * 3 * MAX_H);
		pv.g_lastX = pv.g_lastM + MAX_H;
		pv.g_lastY = pv.g_lastX + MAX_H;
		pv.chunk = mem.chunk;
	}
	__syncthreads();
	//////////////////////////////////////////////////////////////////////

	while (true)
	{
		if (TID == 0)
			compIndex = atomicAdd(g_counter, 1);
		__syncthreads();
		if (compIndex >= mem.nres)
			break;
		__syncthreads();



		if (TID == 0)
		{
			pv.rs = mem.r[mem.cmpR[compIndex]];
			pv.hap = mem.h[mem.cmpH[compIndex]];
			pv.buffsz = 0;
			pv.buffstart = 0;
			pv.H = pv.hap.H;
			pv.R = pv.rs.R;
			pv.nblockrows = ((pv.R+1)+(ROWS-1)) / ROWS;
			pv.nblockcols = (((ROWS)+(pv.H+1)-1) + (DIAGS-1)) / DIAGS;
			pv.m = pv.m1; pv.mp = pv.m2; pv.mpp = pv.m3;
			pv.y = pv.y1; pv.yp = pv.y2; pv.ypp = pv.y3;
			pv.x = pv.x1; pv.xp = pv.x2; pv.xpp = pv.x3;
		}
		__syncthreads();

		for (i=0; i < pv.nblockrows; i++)
		{
			load_read_data(i, &pv);
			__syncthreads();

			if (TID == 0)
			{
				pv.buffstart=0;
				pv.buffsz=0;
			}
			__syncthreads();

			for (j=0; j < pv.nblockcols; j++)
			{
				block<ROWS, DIAGS, BUFFSIZE>(i, j, &pv);
				__syncthreads();
			}
			__syncthreads();

			for (k=TID; k < pv.buffsz; k+=NTH)
			{
				pv.g_lastM[pv.buffstart + k] = pv.buffM[k];
				pv.g_lastX[pv.buffstart + k] = pv.buffX[k];
				pv.g_lastY[pv.buffstart + k] = pv.buffY[k];
			}
			__syncthreads();
		}

		if (TID == 0)
		{
			int resr = (int)pv.R - (int)(pv.nblockrows-1) * (int)ROWS;
			int swaps = ((int)DIAGS * (int)pv.nblockcols -resr - (int)(pv.H + 1)) % 3;
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

void init_cuda_memory(Memory *h_mem, Memory *d_mem)
{
	*d_mem = *h_mem;

	d_mem->g = NULL;
	d_mem->phred_to_prob = NULL;

	cudaError err;

	cudaMalloc((void **) &(d_mem->r), d_mem->nr * sizeof(ReadSequence));
	err = cudaThreadSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(0);
	}
	cudaMemcpy(d_mem->r, h_mem->r, d_mem->nr * sizeof(ReadSequence), cudaMemcpyHostToDevice);

	cudaMalloc((void **) &(d_mem->h), d_mem->nh * sizeof(Haplotype));
	err = cudaThreadSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(0);
	}
	cudaMemcpy(d_mem->h, h_mem->h, d_mem->nh * sizeof(Haplotype), cudaMemcpyHostToDevice);

	cudaMalloc((void **) &(d_mem->chunk), d_mem->chunk_sz);
	err = cudaThreadSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(0);
	}
	cudaMemcpy(d_mem->chunk, h_mem->chunk, d_mem->chunk_sz, cudaMemcpyHostToDevice);

	cudaMalloc((void **) &(d_mem->cmpH), h_mem->nres * sizeof(unsigned long));
	err = cudaThreadSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(0);
	}
	cudaMemcpy(d_mem->cmpH, h_mem->cmpH, d_mem->nres * sizeof(unsigned long), cudaMemcpyHostToDevice);

	cudaMalloc((void **) &(d_mem->cmpR), h_mem->nres * sizeof(unsigned long));
	err = cudaThreadSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(0);
	}
	cudaMemcpy(d_mem->cmpR, h_mem->cmpR, h_mem->nres * sizeof(unsigned long), cudaMemcpyHostToDevice);

	cudaMalloc((void **) &(d_mem->res), h_mem->nres * sizeof(double));
	err = cudaThreadSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));
		exit(0);
	}
}

void cuda_free_memory(Memory m)
{
	cudaFree(m.h);
	cudaFree(m.r);
	cudaFree(m.chunk);
	cudaFree(m.cmpH);
	cudaFree(m.cmpR);
	cudaFree(m.res);
}

int split_allocate(Memory *big, Memory *small)
{
	static ul lastGroup = 0;
	static ul offset_res = 0;
	ul offset_h = big->g[lastGroup].fstH;
	ul offset_r = big->g[lastGroup].fstR;
	ul chunk_begin, chunk_end;
	ul j;
	ul fstG = lastGroup;

	if (lastGroup >= big->ng)
		return 1;

	small->nres=0;
	small->ng=0;
	small->nh=0;
	small->nr=0;
	small->chunk_sz=0;
	while ((lastGroup < big->ng) && (small->nres + big->g[lastGroup].nR * big->g[lastGroup].nH < MAX_COMPARISONS_PER_SPLIT))
	{
		small->nres += big->g[lastGroup].nR * big->g[lastGroup].nH;
		small->ng++;
		small->nh += big->g[lastGroup].nH;
		small->nr += big->g[lastGroup].nR;
		lastGroup++;
	}
	if (small->nres == 0)
	{
		fprintf(stderr, "There exists a group with more than MAX_COMPARISONS_PER_SPLIT comparisons\n");
		exit(0);
	}

	chunk_begin = big->r[big->g[fstG].fstR].r;
	chunk_end = big->h[big->g[lastGroup-1].fstH + big->g[lastGroup-1].nH - 1].h + big->h[big->g[lastGroup-1].fstH + big->g[lastGroup-1].nH - 1].H +1;
	small->chunk_sz = (chunk_end - chunk_begin + 1) + 16;
	small->h = (Haplotype *) malloc(small->nh * sizeof(Haplotype));
	small->r = (ReadSequence *) malloc(small->nr * sizeof(ReadSequence));
	small->cmpH = (unsigned long *) malloc(small->nres * sizeof(unsigned long));
	small->cmpR = (unsigned long *) malloc(small->nres * sizeof(unsigned long));
	small->chunk = (char *) malloc(small->chunk_sz);

	memcpy(small->chunk, big->chunk + chunk_begin, chunk_end - chunk_begin + 1);

	for (j = 0; j < small->nh; j++)
	{
		small->h[j].H = big->h[offset_h + j].H;
		small->h[j].h = big->h[offset_h + j].h - chunk_begin;
	}
	for (j = 0; j < small->nr; j++)
	{
		small->r[j].R = big->r[offset_r + j].R;
		small->r[j].r = big->r[offset_r + j].r - chunk_begin;
		small->r[j].qual = big->r[offset_r + j].qual - chunk_begin;
		small->r[j].ins = big->r[offset_r + j].ins - chunk_begin;
		small->r[j].del = big->r[offset_r + j].del - chunk_begin;
		small->r[j].cont = big->r[offset_r + j].cont - chunk_begin;
	}
	for (j = 0; j < small->nres; j++)
	{
		small->cmpH[j] = big->cmpH[offset_res + j] - offset_h;
		small->cmpR[j] = big->cmpR[offset_res + j] - offset_r;
	}
	offset_res += small->nres;
	return 0;
}

void split_free(Memory m)
{
	free(m.h);
	free(m.r);
	free(m.chunk);
	free(m.cmpH);
	free(m.cmpR);
}

void showtimes(Timers &t)
{
	printf("INIT_MEMORY_TIME %e\n", 1000.0 * (t.t0 - t.start));
	printf("CUDA_SUBTIME_SPLIT %e\n", 1000.0 * t.split);
	printf("CUDA_SUBTIME_COPYING_H_TO_D %e\n", 1000.0 * t.HtoD);
	printf("CUDA_SUBTIME_COMPUTATION %e\n", 1000.0 * t.compute);
	printf("CUDA_SUBTIME_COPYING_D_TO_H %e\n", 1000.0 * t.DtoH);
	printf("COMPUTATION_TIME %e\n", 1000.0 * (t.t6 - t.t0));
	printf("OUTPUT_TIME %e\n", 1000.0 * (t.t7 - t.t6));
	printf("TOTAL_TIME %e\n", 1000.0 * (t.t7 - t.start));
}

void check_args(int argc, char **argv)
{
	(void) argv;
	if (argc != 3)
	{
		printf("\nUsage: <binary> <input file> <output>\n\n");
		exit(0);
	}
}

int main(int argc, char **argv)
{
	Memory hbig, hsmall, dev;
	ul already = 0;
	dbl *g_lastlines;
	int *g_counter, counter = 0;
	Timers t;

	t.start = right_now();
	init_memory(argv[1], &hbig);
	cudaMalloc(&g_lastlines, 3 * sizeof(dbl) * MAX_H * BLOCKS_PER_SM * SMS * 2);
	cudaMalloc(&g_counter, sizeof(int));

	dim3 gridDim(BLOCKS_PER_SM * SMS);
	dim3 blockDim(32, 3);

	t.t0 = right_now();
	while (0 == 0)
	{
		t.t1 = right_now();
		if (split_allocate(&hbig, &hsmall) != 0)
			break;
		t.t2 = right_now();
		cudaThreadSynchronize();
		init_cuda_memory(&hsmall, &dev);
		t.t3 = right_now();
		cudaMemcpy(g_counter, &counter, sizeof(int), cudaMemcpyHostToDevice);
		compare<32, 60, 30><<<gridDim, blockDim>>>(dev, g_lastlines, g_counter);
		cudaThreadSynchronize();
		t.t4 = right_now();
		cudaMemcpy(hbig.res + already, dev.res, dev.nres * sizeof(dbl), cudaMemcpyDeviceToHost);
		cudaThreadSynchronize();
		t.t5 = right_now();
		already += dev.nres;
		cuda_free_memory(dev);
		split_free(hsmall);
		cudaThreadSynchronize();
		t.acc();
	}
	cudaFree(g_lastlines);
	cudaFree(g_counter);

	t.t6 = right_now();
	output(hbig.res, hbig.nres, argv[2]);
	t.t7 = right_now();

	showtimes(t);
	return 0;
}

