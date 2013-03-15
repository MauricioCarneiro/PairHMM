#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "common.h"

/*
	IMPORTANT: IF THERE EXISTS A GROUP WITH MORE COMPARISONS THAN THE
	MAX_COMPARISONS_PER_SPLIT BOUND, THEN IT'S GOING TO HAPPEN A SEGMENTATION FAULT
*/
#define MAX_COMPARISONS_PER_SPLIT ((unsigned long) 50000)

#define MAX_R 256
#define MAX_H 640

#define PRECOMPUTED_CONSTANTS

typedef unsigned long ul;
typedef unsigned char uc;
typedef double dbl;

__device__ inline dbl __m(dbl mp, dbl xp, dbl yp, dbl MM, dbl GapM, uc r, uc h, dbl Q)
{
	return ((r==h || r=='N' || h=='N') ? 1.0-Q : Q) * (mp*MM+xp*GapM+yp*GapM);
}

__device__ inline dbl __x(dbl mp, dbl xp, dbl MX, dbl XX)
{
	return mp * MX + xp * XX;
}

__device__ inline dbl __y(dbl mp, dbl yp, dbl MY, dbl YY)
{
	return mp * MY + yp * YY;
}

__global__ void go(Memory mem)
{
	ul tx = threadIdx.x;
	ul ty = threadIdx.y;
	ul tt = tx + blockDim.x * ty;
	ul nthx = blockDim.x;
	ul nth = blockDim.x * blockDim.y;
	ul compIndex = (ul) (blockIdx.x + gridDim.x * blockIdx.y);

	if (compIndex >= mem.nres)
		return;

	#ifdef PRECOMPUTED_CONSTANTS
	__shared__ dbl MM[MAX_R];
	__shared__ dbl MX[MAX_R];
	__shared__ dbl MY[MAX_R];
	__shared__ dbl GapM[MAX_R];
	__shared__ dbl XX[MAX_R];
	__shared__ dbl YY[MAX_R];
	__shared__ dbl Qual[MAX_R];
	#else
    /*dbl MM, MX, MY, GapM, XX, YY, Qual;*/
	#endif
	__shared__ dbl m1[MAX_R], m2[MAX_R], m3[MAX_R];
	__shared__ dbl x1[MAX_R], x2[MAX_R], x3[MAX_R];
	__shared__ dbl y1[MAX_R], y2[MAX_R], y3[MAX_R];
	__shared__ dbl ph2pr[128];
	__shared__ char _h[MAX_H], _r[MAX_R];
	__shared__ uc _q[MAX_R], _i[MAX_R], _d[MAX_R], _c[MAX_R];

	ul r, c, d, rfrom, rto, cfrom, cto;
	ul i, R, H, cmph, cmpr;
	dbl *sw;
	dbl *m, *mp, *mpp, *y, *yp, *ypp, *x, *xp, *xpp;
	ReadSequence rs;
	Haplotype hap;

	m = m1; mp = m2; mpp = m3;
	x = x1; xp = x2; xpp = x3;
	y = y1; yp = y2; ypp = y3;
	cmph = mem.cmpH[compIndex];
	cmpr = mem.cmpR[compIndex];
	rs = mem.r[cmpr];
	hap = mem.h[cmph];
	H = hap.H;
	R = rs.R;

	for (i = tt; i < R; i+=nth)
	{
		_r[i] = mem.chunk[rs.r + i];
		_q[i] = mem.chunk[rs.qual + i];
		_i[i] = mem.chunk[rs.ins + i];
		_d[i] = mem.chunk[rs.del + i];
		_c[i] = mem.chunk[rs.cont + i];
	}

	for (i = tt; i < H; i+=nth)
		_h[i] = mem.chunk[hap.h + i];

	for (i = tt; i < 128; i+=nth)
		ph2pr[i] = pow(10.0, -((double) i) / 10.0);

	__syncthreads();

	#ifdef PRECOMPUTED_CONSTANTS
	for (i = tt; i < R; i+=nth)
	{
		MM[i+1] = 1.0 - ph2pr[(_i[i] + _d[i]) & 127];
		MX[i+1] = ph2pr[_i[i]];
		MY[i+1] = ph2pr[_d[i]];
		GapM[i+1] = 1.0 - ph2pr[_c[i]];
		XX[i+1] = ph2pr[_c[i]];
		YY[i+1] = ph2pr[_c[i]];
		Qual[i+1] = ph2pr[_q[i]];
	}
	if (tt == 0)
	{
		MM[0] = 1.0 - ph2pr[90];
		MX[0] = ph2pr[45];
		MY[0] = 1.0;
		GapM[0] = 1.0 - ph2pr[10];
		XX[0] = ph2pr[10];
		YY[0] = 1.0;
		Qual[0] = 0.0; //does not matter...
		MY[R] = 1.0;
		YY[R] = 1.0;
	}
	#endif


	// d = 0
	if (tt == 0)
	{
		/*m[0] = 1.0;*/
		m[0] = 1.0e300 / (H > R ? H - R + 1 : 1);
		x[0] = 0.0;
		y[0] = 0.0;
	}

	for (d=1; d < R+H+1; d++)
	{
		__syncthreads();

		sw=mpp; mpp=mp; mp=m; m=sw; 
		sw=xpp; xpp=xp; xp=x; x=sw;
		sw=ypp; ypp=yp; yp=y; y=sw;

		rfrom = (d < H) ? 0 : d - H;
		rto = (d > R) ? R : d;
		cfrom = (d < H) ? d : H;
		cto = (d > R) ? d - R : 0;

		if (rfrom == 0) 
		{
			if (tt == 0)
			{
				m[0] = 0.0;
				x[0] = 0.0;
				y[0] = mp[0] + yp[0];
			}
			rfrom++;
			cfrom--;
		}

		if (cto == 0)
		{
			if (tt == 0)
			{
				m[d] = 0.0;
				x[d] = mp[rto-1]*ph2pr[_i[rto-1]] + xp[rto-1]*ph2pr[_c[rto-1]];
				y[d] = 0.0;
			}
			cto++;
			rto--;
		}

		for (r=rfrom+tx, c=cfrom-tx; r <= rto && c >= cto; r+=nthx, c-=nthx)
			if (ty == 0)
			{
				#ifdef PRECOMPUTED_CONSTANTS
				m[r] = __m(mpp[r-1], xpp[r-1], ypp[r-1], MM[r], GapM[r], _r[r-1], _h[c-1], Qual[r]);
				#else
				//MM = 1.0 - ph2pr[(_i[r-1] + _d[r-1]) & 127];
				//GapM = 1.0 - ph2pr[_c[r-1]];
				//Qual = ph2pr[_q[r-1]];
				m[r] = __m(mpp[r-1], xpp[r-1], ypp[r-1], 1.0 - ph2pr[(_i[r-1] + _d[r-1]) & 127], 1.0 - ph2pr[_c[r-1]], _r[r-1], _h[c-1], ph2pr[_q[r-1]]);
				#endif
			}
			else if (ty == 1)
			{
				#ifdef PRECOMPUTED_CONSTANTS
				x[r] = __x(mp[r-1], xp[r-1], MX[r], XX[r]);
				#else
				//MX = ph2pr[_i[r-1]];
				//XX = ph2pr[_c[r-1]];
				x[r] = __x(mp[r-1], xp[r-1], ph2pr[_i[r-1]], ph2pr[_c[r-1]]);
				#endif
			}
			else if (ty == 2)
			{
				#ifdef PRECOMPUTED_CONSTANTS
				y[r] = __y(mp[r], yp[r], MY[r], YY[r]);
				#else
				//MY = (r == R) ? 1.0 : ph2pr[_d[r-1]];
				//YY = (r == R) ? 1.0 : ph2pr[_c[r-1]];
				y[r] = (r == R) ? __y(mp[r], yp[r], 1.0, 1.0) : __y(mp[r], yp[r], ph2pr[_d[r-1]], ph2pr[_c[r-1]]);
				#endif
			}
	}
	__syncthreads();

	if (tt == 0)
		mem.res[compIndex] = log10(m[R] + x[R] + y[R]) - 300.0;

	return;
}

void init_cuda_memory(Memory *h_mem, Memory *d_mem)
{
	*d_mem = *h_mem;

	d_mem->g = NULL;
	d_mem->phred_to_prob = NULL;

	if (cudaMalloc((void **) &(d_mem->r), d_mem->nr * sizeof(ReadSequence)) != 0) { printf("Error in cuda memory allocation, line %d\n", __LINE__); exit(0); }
	cudaMemcpy(d_mem->r, h_mem->r, d_mem->nr * sizeof(ReadSequence), cudaMemcpyHostToDevice);
	if (cudaMalloc((void **) &(d_mem->h), d_mem->nh * sizeof(Haplotype)) != 0) { printf("Error in cuda memory allocation, line %d\n", __LINE__); exit(0); }
	cudaMemcpy(d_mem->h, h_mem->h, d_mem->nh * sizeof(Haplotype), cudaMemcpyHostToDevice);
	if (cudaMalloc((void **) &(d_mem->chunk), d_mem->chunk_sz) != 0) { printf("Error in cuda memory allocation, line %d\n", __LINE__); exit(0); }
	cudaMemcpy(d_mem->chunk, h_mem->chunk, d_mem->chunk_sz, cudaMemcpyHostToDevice);
	if (cudaMalloc((void **) &(d_mem->cmpH), h_mem->nres * sizeof(unsigned long)) != 0) { printf("Error in cuda memory allocation, line %d\n", __LINE__); exit(0); }
	cudaMemcpy(d_mem->cmpH, h_mem->cmpH, d_mem->nres * sizeof(unsigned long), cudaMemcpyHostToDevice);
	if (cudaMalloc((void **) &(d_mem->cmpR), h_mem->nres * sizeof(unsigned long)) != 0) { printf("Error in cuda memory allocation, line %d\n", __LINE__); exit(0); }
	cudaMemcpy(d_mem->cmpR, h_mem->cmpR, h_mem->nres * sizeof(unsigned long), cudaMemcpyHostToDevice);
	if (cudaMalloc((void **) &(d_mem->res), h_mem->nres * sizeof(double)) != 0) { printf("Error in cuda memory allocation, line %d\n", __LINE__); exit(0); }
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
	static unsigned long lastGroup = 0;
	static unsigned long offset_res = 0;
	unsigned long offset_h = big->g[lastGroup].fstH;
	unsigned long offset_r = big->g[lastGroup].fstR;
	unsigned long chunk_begin, chunk_end;
	unsigned long j;
	unsigned long fstG = lastGroup;

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

/*
void show_memory(Memory m)
{
	unsigned long x;
	printf("Showing memory...\n");
	printf("m.chunk = ***\n");
	for (x = 0; x < m.chunk_sz; x++)
		printf("%3d ", x);
	printf("\n");
	for (x = 0; x < m.chunk_sz; x++)
		printf("%3d ", (int) m.chunk[x]);
	printf("\n***\n");
}
*/

int main(int argc, char **argv)
{
	Memory hbig, hsmall, dev;
	double start, t0, t1, t2, t3, t4, t5, t6, t7;
	double splitting = 0, computing = 0, copyingHtoD = 0, copyingDtoH = 0;
	unsigned long already = 0;
	if (argc != 3)
	{
		printf("\nUsage: <binary> <input file> <output>\n\n");
		exit(0);
	}

	start = right_now();
	init_memory(argv[1], &hbig);
	t0 = right_now();
	while (0 == 0)
	{
		t1 = right_now();
		if (split_allocate(&hbig, &hsmall) != 0) break;
		t2 = right_now();
		init_cuda_memory(&hsmall, &dev);
		t3 = right_now();
		dim3 gridDim(65000, (dev.nres + 65000 - 1) / 65000);
		dim3 blockDim(128,3);
		go<<<gridDim, blockDim>>>(dev);
		cudaThreadSynchronize();
		t4 = right_now();
		cudaMemcpy(hbig.res + already, dev.res, dev.nres * sizeof(double), cudaMemcpyDeviceToHost);
		cudaThreadSynchronize();
		t5 = right_now();
		already += dev.nres;
		cuda_free_memory(dev);
		split_free(hsmall);
		cudaThreadSynchronize();

		splitting += (t2 - t1);
		copyingHtoD += (t3 - t2);
		computing += (t4 - t3);
		copyingDtoH += (t5 - t4);
	}
	t6 = right_now();
	output(hbig.res, hbig.nres, argv[2]);
	t7 = right_now();
	printf("INIT_MEMORY_TIME %e\n", 1000.0 * (t0 - start));
	printf("CUDA_SUBTIME_SPLIT %e\n", 1000.0 * splitting);
	printf("CUDA_SUBTIME_COPYING_H_TO_D %e\n", 1000.0 * copyingHtoD);
	printf("CUDA_SUBTIME_COMPUTATION %e\n", 1000.0 * computing);
	printf("CUDA_SUBTIME_COPYING_D_TO_H %e\n", 1000.0 * copyingDtoH);
	printf("COMPUTATION_TIME %e\n", 1000.0 * (t6 - t0));
	printf("OUTPUT_TIME %e\n", 1000.0 * (t7 - t6));
	printf("TOTAL_TIME %e\n", 1000.0 * (t7 - start));
	return 0;
}

