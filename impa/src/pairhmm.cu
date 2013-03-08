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

__global__ void go2(Memory mem)
{
	__shared__ double p[MAX_R][6];
	__shared__ double m1[MAX_R], m2[MAX_R], m3[MAX_R];
	__shared__ double x1[MAX_R], x2[MAX_R], x3[MAX_R];
	__shared__ double y1[MAX_R], y2[MAX_R], y3[MAX_R];
	__shared__ double ph2pr[256];
	__shared__ char h[MAX_H];
	__shared__ char r[MAX_R], qual[MAX_R], ins[MAX_R], del[MAX_R], cont[MAX_R];

    unsigned long tx = threadIdx.x;
    unsigned long i = (unsigned long) blockIdx.x + (unsigned long) 65000 * (unsigned long) blockIdx.y;
	unsigned long row, col, diag, i1, R, H, hIndex, rIndex;
	double *m = m1, *mp = m2, *mpp = m3;
	double *y = y1, *yp = y2, *ypp = y3;
	double *x = x1, *xp = x2, *xpp = x3;
	double *swapper;

	if (i >= mem.nres) 
		return;

    hIndex = mem.cmpH[i];
    rIndex = mem.cmpR[i];
    H = mem.h[hIndex].H;
    R = mem.r[rIndex].R;

    /* Collaborative Memory Inicialization -- BEGIN */

    /* read sequence, qual, ins, del, cont */
    for (i1 = tx; i1 < R; i1+=32)
    {
        r[i1] = mem.chunk[mem.r[rIndex].r + i1];
        qual[i1] = mem.chunk[mem.r[rIndex].qual + i1];
        ins[i1] = mem.chunk[mem.r[rIndex].ins + i1];
        del[i1] = mem.chunk[mem.r[rIndex].del + i1];
        cont[i1] = mem.chunk[mem.r[rIndex].cont + i1];
    }

    /* haplotype */
    for (i1 = tx; i1 < H; i1+=32)
        h[i1] = mem.chunk[mem.h[hIndex].h + i1];

    /* ph2pr */
    for (i1 = tx; i1 < 256; i1+=32)
        ph2pr[i1] = pow(10.0, -((double) i1) / 10.0);

    /* transitions */
    for (i1 = tx; i1 < R+1; i1+=32)
        if (i1 == 0)
        {
            p[i1][MtoM] = 1.0 - ph2pr[90];
            p[i1][GapToM] = 1.0 - ph2pr[10];
            p[i1][MtoX] = ph2pr[45];
            p[i1][XtoX] = ph2pr[10];
            p[i1][MtoY] = 1.0;
            p[i1][YtoY] = 1.0;
        }
        else
        {
            p[i1][MtoM] = 1.0 - ph2pr[(ins[i1-1] & 127) + (del[i1-1] & 127)];
            p[i1][GapToM] = 1.0 - ph2pr[cont[i1-1] & 127];
            p[i1][MtoX] = ph2pr[ins[i1-1] & 127];
            p[i1][XtoX] = ph2pr[cont[i1-1] & 127];
            p[i1][MtoY] = (i1 == R) ? 1.0 : ph2pr[del[i1-1] & 127];
            p[i1][YtoY] = (i1 == R) ? 1.0 : ph2pr[cont[i1-1] & 127];
        }

    /* Collaborative Memory Inicialization -- END */

    /* Collaborative Computation -- BEGIN */

    /* diagonal 0 */
    /*m[0] = 1.0;*/
	m[0] = exp10(300.0) / (H > R ? H - R + 1 : 1); 
    x[0] = 0.0;
    y[0] = 0.0;

    /* diagonal diag */
    for (diag = 1; diag < R+H+1; diag++)
    {
        swapper = mpp; mpp = mp; mp = m; m = swapper;
        swapper = xpp; xpp = xp; xp = x; x = swapper;
        swapper = ypp; ypp = yp; yp = y; y = swapper;

        for (row = tx; row <= R; row+=32)
            if ((diag <= H + row) && (diag >= row))
            {
				col = diag - row;
                m[row] = ((row == 0) || (col == 0)) ? 0.0 : ( (r[row-1] == h[col-1] || r[row-1] == 'N' || h[col-1] == 'N') ? 1.0 - ph2pr[qual[row-1] & 127] : ph2pr[qual[row-1] & 127] ) * (mpp[row-1] * p[row][MtoM] + xpp[row-1] * p[row][GapToM] + ypp[row-1] * p[row][GapToM]);
                x[row] = (row == 0) ? 0.0 : mp[row-1] * p[row][MtoX] + xp[row-1] * p[row][XtoX];
                y[row] = (col == 0) ? 0.0 : mp[row] * p[row][MtoY] + yp[row] * p[row][YtoY];
            }
    }

    /* Collaborative Computation -- END */

	if (tx == 0)
	{
    	/*mem.res[i] = log10(m[R] + x[R] + y[R]);*/
    	mem.res[i] = log10(m[R] + x[R] + y[R]) - 300.0;
	}

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
		go2<<<gridDim, 32>>>(dev);
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

