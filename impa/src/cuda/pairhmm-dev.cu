#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <ctime>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>

#define BIGGEST_NUM double

#define COMPDIAGS 60
#define COMPBUFFSIZE 30

#define MIN(a,b) (((a)<(b))?(a):(b))
#define TX (threadIdx.x)
#define TY (threadIdx.y)
#define TID (threadIdx.x + BLOCKWIDTH * threadIdx.y)
#define NTH (BLOCKHEIGHT * BLOCKWIDTH)

#define MAX_COMPARISONS_PER_SPLIT ((unsigned long) 35000)
#define MAX_H 5120
#define MAX_R 1536
#define MAX_SIMULTANEOUS_BLOCKS 150

typedef struct
{
	unsigned long nR, nH, fstR, fstH;
} Group;

typedef struct
{
	unsigned long R;
	int r, qidc;
} ReadSequence;

typedef struct
{
	unsigned long H;
	int h;
} Haplotype;

typedef struct
{
	Group *g;
	Haplotype *h;
	ReadSequence *r;
	char *chunk;
	unsigned long *cmpH;
	unsigned long *cmpR;
	BIGGEST_NUM *res;
	char *flag;
	unsigned long ng, nh, nr, chunk_sz, nres;
} Memory;

typedef unsigned long ul;
typedef unsigned char uc;

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
struct PubVars {
	ul H, R;
	char h[ROWS+DIAGS-1];
	char r[ROWS];
	char4 qidc[ROWS];
	ul nblockrows, nblockcols;
	Haplotype hap;
	ReadSequence rs;
	char *chunk;

	NUM ph2pr[128];
	NUM m1[ROWS], m2[ROWS], m3[ROWS];
	NUM x1[ROWS], x2[ROWS], x3[ROWS];
	NUM y1[ROWS], y2[ROWS], y3[ROWS];
	NUM *m, *mp, *mpp;
	NUM *x, *xp, *xpp;
	NUM *y, *yp, *ypp;
	NUM *g_lastM, *g_lastX, *g_lastY;
	NUM lastM[DIAGS+1], lastX[DIAGS+1], lastY[DIAGS+1];
	NUM buffM[BUFFSIZE], buffX[BUFFSIZE], buffY[BUFFSIZE];
	ul buffsz;
	ul buffstart;
	NUM result;
};

double right_now(void);
void init_memory(const char *fn, Memory *m);
void output(BIGGEST_NUM *r, unsigned long nr, const char *filename);

template<class T>
__device__ static inline T INITIAL_CONSTANT();

template<>
__device__ static inline float INITIAL_CONSTANT<float>() {
	return 1e32f;
}

template<>
__device__ static inline double INITIAL_CONSTANT<double>() {
	return ldexp(1.0, 1020);
}

template<class T>
__device__ static inline T MIN_ACCEPTED();

template<>
__device__ static inline float MIN_ACCEPTED<float>() {
	return 1e-28f;
}

template<>
__device__ static inline double MIN_ACCEPTED<double>() {
	return 0.0;
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void flush_buffer( PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	for (int k = TID; k < pv.buffsz; k += NTH) {
		pv.g_lastM[pv.buffstart + k] = pv.buffM[k];
		pv.g_lastX[pv.buffstart + k] = pv.buffX[k];
		pv.g_lastY[pv.buffstart + k] = pv.buffY[k];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void load_previous_results( int &j, PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	for (int c = ((j == 0) ? 1 : 0) + TID; c < MIN(pv.H + 1 - j * DIAGS, DIAGS) + 1; c += NTH) {
		pv.lastM[c] = pv.g_lastM[j * DIAGS + c - 1];
		pv.lastX[c] = pv.g_lastX[j * DIAGS + c - 1];
		pv.lastY[c] = pv.g_lastY[j * DIAGS + c - 1];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void load_hap_data( int &j, PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	if (TID < DIAGS)
		for (int c = TID; c < ROWS - 1; c += DIAGS)
			pv.h[c] = pv.h[c + DIAGS];

	__syncthreads();

	for (int c = ((j == 0) ? 1 : 0) + TID; c < MIN(DIAGS, pv.H + 1 - j * DIAGS); c += NTH)
		pv.h[ROWS - 1 + c] = pv.chunk[pv.hap.h + j * DIAGS + c - 1];
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void load_read_data( int &i, PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	for (int r = ((i == 0) ? 1 : 0) + TID; r < MIN(ROWS, pv.R + 1 - i * ROWS); r += NTH) {
		pv.r[r] = pv.chunk[pv.rs.r + i * ROWS + r - 1];
		pv.qidc[r].x = pv.chunk[pv.rs.qidc + 4 * (i * ROWS + r - 1) + 0];
		pv.qidc[r].y = pv.chunk[pv.rs.qidc + 4 * (i * ROWS + r - 1) + 1];
		pv.qidc[r].z = pv.chunk[pv.rs.qidc + 4 * (i * ROWS + r - 1) + 2];
		pv.qidc[r].w = pv.chunk[pv.rs.qidc + 4 * (i * ROWS + r - 1) + 3];
	}
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void notfirstline_firstcolum( int &r, PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv, NUM &_m, NUM &_x ) {
	pv.m[r] = (NUM(0.0));
	pv.x[r] = _m * pv.ph2pr[pv.qidc[r].y] + _x * pv.ph2pr[pv.qidc[r].w];
	pv.y[r] = (NUM(0.0));
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void notfirstline_notfirstcolumn( int &r, int &i, int &diag, PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv, NUM &M_m, NUM &M_x, NUM &M_y, NUM &X_m, NUM &X_x, NUM &Y_m, NUM &Y_y ) {
	NUM t1, t2, t3, dist;

	t1 = (NUM(1.0)) - pv.ph2pr[(pv.qidc[r].y + pv.qidc[r].z) & 127];
	t2 = (NUM(1.0)) - pv.ph2pr[pv.qidc[r].w];
	t3 = pv.ph2pr[pv.qidc[r].x];

	dist = (pv.r[r] == pv.h[ROWS - 1 + diag - r] || pv.r[r] == 'N' || pv.h[ROWS - 1 + diag - r] == 'N') ? (NUM(1.0)) - t3 : t3;
		pv.m[r] = dist * (M_m * t1 + M_x * t2 + M_y * t2);

	t1 = pv.ph2pr[pv.qidc[r].y];
	t2 = pv.ph2pr[pv.qidc[r].w];
	pv.x[r] = X_m * t1 + X_x * t2;

	t1 = ((unsigned)r + i * ROWS == pv.R) ? (NUM(1.0)) : pv.ph2pr[pv.qidc[r].z];
	t2 = ((unsigned)r + i * ROWS == pv.R) ? (NUM(1.0)) : pv.ph2pr[pv.qidc[r].w];
	pv.y[r] = Y_m * t1 + Y_y * t2;
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void firstline_firstcolum( PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	pv.m[0] = (NUM(0.0));
	pv.x[0] = (NUM(0.0));
	pv.y[0] = INITIAL_CONSTANT<NUM>() / pv.H;
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void firstline_notfirstcolum( PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	pv.m[0] = (NUM(0.0));
	pv.x[0] = (NUM(0.0));
	pv.y[0] = INITIAL_CONSTANT<NUM>() / pv.H;
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void rotatediags(	PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	NUM *sw;
	sw = pv.mpp; pv.mpp = pv.mp; pv.mp = pv.m; pv.m = sw;
	sw = pv.xpp; pv.xpp = pv.xp; pv.xp = pv.x; pv.x = sw;
	sw = pv.ypp; pv.ypp = pv.yp; pv.yp = pv.y; pv.y = sw;
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__device__ inline void block( int &i, int &j, PubVars<ROWS, DIAGS, BUFFSIZE, NUM> &pv ) {
	if (i > 0)
		load_previous_results(j, pv);

	load_hap_data(j, pv);

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
							firstline_firstcolum(pv);
						else
							firstline_notfirstcolum(pv);
					else
						if (c == 0)
							notfirstline_firstcolum(r, pv, pv.lastM[1], pv.lastX[1]);
						else
							notfirstline_notfirstcolumn(r, i, diag, pv, pv.lastM[diag], pv.lastX[diag], pv.lastY[diag], pv.lastM[diag + 1], pv.lastX[diag + 1], pv.mp[r], pv.yp[r]);
				else
					if (c == 0)
						notfirstline_firstcolum(r, pv, pv.mp[r-1], pv.xp[r-1]);
					else
						notfirstline_notfirstcolumn(r, i, diag, pv, pv.mpp[r-1], pv.xpp[r-1], pv.ypp[r-1], pv.mp[r-1], pv.xp[r-1], pv.mp[r], pv.yp[r]);
		}

		__syncthreads();

		r = nRows - 1;
		c = diag - r + j * DIAGS;

		if ((TID == 0) && (c >= 0) && (c <= (int)pv.H))
			if (i < pv.nblockrows - 1)
			{
				pv.buffM[pv.buffsz] = pv.m[r];
				pv.buffX[pv.buffsz] = pv.x[r];
				pv.buffY[pv.buffsz] = pv.y[r];
				pv.buffsz++;
			}
			else
				pv.result += (pv.m[r] + pv.x[r]);

		__syncthreads();

		if (pv.buffsz == BUFFSIZE)
		{
			flush_buffer(pv);

			if (TID == 0)
			{
				pv.buffstart += BUFFSIZE;
				pv.buffsz = 0;
			}
			__syncthreads();
		}

		if (TID == 0)
			rotatediags(pv);

		__syncthreads();
	}

	return;
}

template <int ROWS, int DIAGS, int BUFFSIZE, typename NUM>
__global__ void compare( Memory mem, NUM *g_lastLinesArr, int *g_lastLinesIndex, int *g_compIndex ) {
	__shared__ PubVars<ROWS, DIAGS, BUFFSIZE, NUM> pv;
	__shared__ int compIndex;
	__shared__ bool flag_zero;

	for (int i = TID; i < 128; i += NTH)
		pv.ph2pr[i] = pow((NUM(10.0)), -(NUM(i)) / (NUM(10.0)));

	if (TID == 0) {
		int lastLinesIndex = atomicAdd(g_lastLinesIndex, 1);
		pv.g_lastM = g_lastLinesArr + (lastLinesIndex * 3 * MAX_H);
		pv.g_lastX = pv.g_lastM + MAX_H;
		pv.g_lastY = pv.g_lastX + MAX_H;
		pv.chunk = mem.chunk;
	}

	__syncthreads();

	for (;;) {
		if (TID == 0) {
			compIndex = atomicAdd(g_compIndex, 1);
			flag_zero = (mem.flag[compIndex] == 0);
		}

		__syncthreads();

		if (compIndex >= mem.nres)
			break;

		if (!flag_zero)
			continue;

		if (TID == 0) {
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
			pv.result = 0;
		}

		__syncthreads();

		for (int i = 0; i < pv.nblockrows; i++) {
			load_read_data(i, pv);

			if (TID == 0) {
				pv.buffstart = 0;
				pv.buffsz = 0;
			}

			__syncthreads();

			for (int j = 0; j < pv.nblockcols; j++) {
				block(i, j, pv);
				__syncthreads();
			}
			flush_buffer(pv);

			__syncthreads();
		}

		if (TID == 0) {
			if (pv.result > MIN_ACCEPTED<NUM>()) {
				mem.flag[compIndex] = 1;
				mem.res[compIndex] = log10(pv.result) - log10(INITIAL_CONSTANT<NUM>());
			}
		}

		__syncthreads();
	}

	return;
}

int split( Memory &h_big, Memory &ret ) {
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

	while ( (lastGroup < h_big.ng) && (ret.nres + h_big.g[lastGroup].nR * h_big.g[lastGroup].nH < MAX_COMPARISONS_PER_SPLIT)) {
		ret.nres += h_big.g[lastGroup].nR * h_big.g[lastGroup].nH;
		ret.ng++;
		ret.nh += h_big.g[lastGroup].nH;
		ret.nr += h_big.g[lastGroup].nR;
		lastGroup++;
	}

	if (ret.nres == 0) {
		fprintf(stderr, "There exists a group with more than MAX_COMPARISONS_PER_SPLIT comparisons\n");
		exit(0);
	}

	chunk_begin = h_big.r[h_big.g[fstG].fstR].r;
	chunk_end = h_big.h[h_big.g[lastGroup-1].fstH + h_big.g[lastGroup-1].nH-1].h + h_big.h[h_big.g[lastGroup-1].fstH + h_big.g[lastGroup-1].nH-1].H + 1;
	ret.chunk_sz = (chunk_end - chunk_begin + 1);
	ret.chunk = h_big.chunk + chunk_begin;

	for (j = 0; j < ret.nh; j++) {
		ret.h[j] = h_big.h[offset_h + j];
		ret.h[j].h -= chunk_begin;
	}

	for (j = 0; j < ret.nr; j++) {
		ret.r[j] = h_big.r[offset_r + j];
		ret.r[j].r -= chunk_begin;
		ret.r[j].qidc -= chunk_begin;
	}

	for (j = 0; j < ret.nres; j++) {
		ret.cmpH[j] = h_big.cmpH[offset_res + j] - offset_h;
		ret.cmpR[j] = h_big.cmpR[offset_res + j] - offset_r;
	}

	offset_res += ret.nres;

	return 0;
}

int main( int argc, char **argv ) {
	Memory h_big, h_small, d_mem;
	ul already = 0;
	void *g_lastlines;
	int *g_compIndex, compIndex = 0;
	int *g_lastLinesIndex, lastLinesIndex = 0;

	struct {
		double start, init_memory, mallocs, comp, output, end, t1, t2;
		float kernel;
	} times;

	times.kernel = 0.f;

	times.start = right_now();
	times.t1 = right_now();
	/*
		init_memory recibe un nombre de archivo y un puntero a una estructura
		Memory, e inicializa h_big (una estructura de tipo Memory) basándose en 
		el contenido del archivo.
	*/
	init_memory(argv[1], &h_big);
	times.t2 = right_now();
	times.init_memory = times.t2 - times.t1;
	times.t1 = right_now();
	h_small.r = (ReadSequence *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(ReadSequence));
	h_small.h = (Haplotype *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(Haplotype));
	h_small.chunk = (char *) malloc(MAX_COMPARISONS_PER_SPLIT * (MAX_H + MAX_R * 5));
	h_small.cmpH = (ul *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	h_small.cmpR = (ul *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	h_small.res = (BIGGEST_NUM *) malloc(MAX_COMPARISONS_PER_SPLIT * sizeof(BIGGEST_NUM));
	h_small.flag = (char *) malloc(MAX_COMPARISONS_PER_SPLIT);
	h_small.g = NULL;

	g_lastlines = NULL;
	g_compIndex = NULL;
	g_lastLinesIndex = NULL;
	d_mem.r = NULL;
	d_mem.h = NULL;
	d_mem.chunk = NULL;
	d_mem.cmpH = NULL;
	d_mem.cmpR = NULL;
	d_mem.res = NULL;

	cudaMalloc(&g_compIndex, sizeof(int));
	cudaMalloc(&g_lastLinesIndex, sizeof(int));
	cudaMalloc(&g_lastlines, 3 * sizeof(BIGGEST_NUM) * MAX_H * MAX_SIMULTANEOUS_BLOCKS);

	cudaMalloc(&(d_mem.r), MAX_COMPARISONS_PER_SPLIT * sizeof(ReadSequence));
	cudaMalloc(&(d_mem.h), MAX_COMPARISONS_PER_SPLIT * sizeof(Haplotype));
	cudaMalloc(&(d_mem.chunk), MAX_COMPARISONS_PER_SPLIT * (MAX_H + MAX_R * 5));
	cudaMalloc(&(d_mem.cmpH), MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	cudaMalloc(&(d_mem.cmpR), MAX_COMPARISONS_PER_SPLIT * sizeof(ul));
	cudaMalloc(&(d_mem.flag), MAX_COMPARISONS_PER_SPLIT);
	cudaMalloc(&(d_mem.res), MAX_COMPARISONS_PER_SPLIT * sizeof(BIGGEST_NUM));
	d_mem.g = NULL;

	if (!g_lastLinesIndex || !g_lastlines || !g_compIndex || !d_mem.r || !d_mem.h || !d_mem.chunk || !d_mem.cmpH || !d_mem.cmpR || !d_mem.res) {
		fprintf(stderr, "Some malloc went wrong...\n");
		exit(0);
	}

	times.t2 = right_now();
	times.mallocs = times.t2 - times.t1;

	times.t1 = right_now();
	while (!split(h_big, h_small)) {
		cudaEvent_t kernel_start, kernel_stop;
		float k_time;
		cudaEventCreate(&kernel_start);
		cudaEventCreate(&kernel_stop);

		d_mem.ng = h_small.ng;
		d_mem.nh = h_small.nh;
		d_mem.nr = h_small.nr;
		d_mem.chunk_sz = h_small.chunk_sz;
		d_mem.nres = h_small.nres;
		cudaMemcpy(d_mem.r, h_small.r, h_small.nr * sizeof(ReadSequence), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.h, h_small.h, h_small.nh * sizeof(Haplotype), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.chunk, h_small.chunk, h_small.chunk_sz, cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.cmpH, h_small.cmpH, h_small.nres * sizeof(ul), cudaMemcpyHostToDevice);
		cudaMemcpy(d_mem.cmpR, h_small.cmpR, h_small.nres * sizeof(ul), cudaMemcpyHostToDevice);

		memset(h_small.flag, 0, MAX_COMPARISONS_PER_SPLIT);
		cudaMemcpy(d_mem.flag, h_small.flag, MAX_COMPARISONS_PER_SPLIT, cudaMemcpyHostToDevice);

		dim3 gridDim(MAX_SIMULTANEOUS_BLOCKS);
		dim3 blockDim(BLOCKWIDTH, BLOCKHEIGHT);
		cudaEventRecord(kernel_start, 0);

		cudaMemcpy(g_compIndex, &compIndex, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(g_lastLinesIndex, &lastLinesIndex, sizeof(int), cudaMemcpyHostToDevice);
		compare<BLOCKWIDTH, COMPDIAGS, COMPBUFFSIZE, float><<<gridDim, blockDim>>>(d_mem, (float *) g_lastlines, g_lastLinesIndex, g_compIndex);

		cudaMemcpy(g_compIndex, &compIndex, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(g_lastLinesIndex, &lastLinesIndex, sizeof(int), cudaMemcpyHostToDevice);
		compare<BLOCKWIDTH, COMPDIAGS, COMPBUFFSIZE, double><<<gridDim, blockDim>>>(d_mem, (double *) g_lastlines, g_lastLinesIndex, g_compIndex);

		cudaEventRecord(kernel_stop, 0);
		cudaEventSynchronize(kernel_stop);

		cudaMemcpy(h_big.res + already, d_mem.res, d_mem.nres * sizeof(BIGGEST_NUM), cudaMemcpyDeviceToHost);
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

double right_now(void) {
	struct timeval v;
	gettimeofday(&v, (struct timezone *) NULL);
	return v.tv_sec + v.tv_usec/1.0e6;
}

void init_memory(const char *fn, Memory *m) {
	char *l = NULL;
	int lastptr;
	size_t sz = 0;
	FILE *i;
	unsigned long nG = 0, nR = 0, nH = 0, tr, th, R = 0, H = 0, x = 0;
	unsigned long y = 0, chunk_sz = 0, nres = 0, curr_g = 0, curr_r = 0;
	unsigned long curr_h = 0, curr_res = 0, gr = 0, r = 0, h = 0;

	std::string qual, ins, del, cont;

	i = fopen(fn, "r");
	while (getline(&l, &sz, i) != -1) {
		if (strcmp(l, "GROUP\n"))
			break;
		nG++;
		if (getline(&l, &sz, i) == -1)
			return;
		sscanf(l, "%lu %lu\n", &tr, &th);
		nR += tr;
		nH += th;
		nres += (tr * th);
		for (x = 0; x < tr; x++) {
			if (getline(&l, &sz, i) == -1) return;
			R = 0; while (l[R] != ' ') R++;
			chunk_sz += (R+1) * 5;
		}
		for (x = 0; x < th; x++) {
			if (getline(&l, &sz, i) == -1) return;
			H = 0; while (l[H] != ' ' && l[H] != '\n') H++;
			chunk_sz += (H+1);
		}
	}

	fclose(i);

	m->ng = nG;
	m->g = (Group *) malloc(m->ng * sizeof(Group));

	m->nh = nH;
	m->h = (Haplotype *) malloc(m->nh * sizeof(Haplotype));

	m->nr = nR;
	m->r = (ReadSequence *) malloc(m->nr * sizeof(ReadSequence));

	m->chunk_sz = chunk_sz;
	m->chunk = (char *) malloc(chunk_sz);
	bzero(m->chunk, chunk_sz);

	m->nres = nres;
	m->res = (BIGGEST_NUM *) malloc(m->nres * sizeof(BIGGEST_NUM));
	m->cmpH = (unsigned long *) malloc(m->nres * sizeof(unsigned long));
	m->cmpR = (unsigned long *) malloc(m->nres * sizeof(unsigned long));

	i = fopen(fn, "r");
	nG = 0;
	curr_g = 0; curr_r = 0; curr_h = 0;
	lastptr = 0;
	while (getline(&l, &sz, i) != -1) {
		if (strcmp(l, "GROUP\n")) break;
		nG++;
		if (getline(&l, &sz, i) == -1) return;

		curr_g = nG - 1;
		sscanf(l, "%lu %lu\n", &(m->g[curr_g].nR), &(m->g[curr_g].nH));

		m->g[curr_g].fstR = curr_r;
		m->g[curr_g].fstH = curr_h;

		for (x = 0; x < m->g[curr_g].nR; x++) {
			if (getline(&l, &sz, i) == -1) 
				return;

			// read the readsequence curr_r

			R = 0; while (l[R] != ' ') R++;
			// R is the number of characters in the read sequence
			m->r[curr_r].R = R; 
			m->r[curr_r].r = lastptr; // m->r[curr_r].r is the position in the chunk in which the readsequence starts. In the chunk, the readsequence has R characters followed by a '\0' character.
			m->r[curr_r].qidc = lastptr + (R+1); // m->r[curr_r].qidc is the position in the chunk in which the readsequence starts
			sscanf(l, "%s ", m->chunk + lastptr); // read r. Put the R characters into the chunk, starting at m->r[curr_r].r
			y = R+1;

			// read qual, ins, del, cont and put them in an interleaved way into the chunk.
			std::istringstream iss(l + y);
			iss >> qual >> ins >> del >> cont;
			std::istringstream issq(qual);
			std::istringstream issi(ins);
			std::istringstream issd(del);
			std::istringstream issc(cont);
			for (int j = 0; j < R; j++) {
				if (issq.peek() == ',') issq.ignore();
				if (issi.peek() == ',') issi.ignore();
				if (issd.peek() == ',') issd.ignore();
				if (issc.peek() == ',') issc.ignore();
				int intq, inti, intd, intc;
				issq >> intq; if (intq < 6) intq = 6;
				issi >> inti;
				issd >> intd;
				issc >> intc;
				m->chunk[lastptr + (R+1) + 4 * j + 0] = (char) intq;
				m->chunk[lastptr + (R+1) + 4 * j + 1] = (char) inti;
				m->chunk[lastptr + (R+1) + 4 * j + 2] = (char) intd;
				m->chunk[lastptr + (R+1) + 4 * j + 3] = (char) intc;
			}

			lastptr += (R+1) * 5;
			curr_r++;
		}
		for (x = 0; x < m->g[curr_g].nH; x++) {
			if (getline(&l, &sz, i) == -1) 
				return;
			H = 0; while (l[H] != ' ' && l[H] != '\n') H++;

			m->h[curr_h].H = H;
			m->h[curr_h].h = lastptr;
			sscanf(l, "%s\n", m->chunk + lastptr);

			lastptr += H+1;
			curr_h++;
		}
	}

	fclose(i);

	curr_res = 0;
	for (gr = 0; gr < m->ng; gr++)
		for (r = 0; r < m->g[gr].nR; r++)
			for (h = 0; h < m->g[gr].nH; h++) {
				m->cmpH[curr_res] = m->g[gr].fstH + h;
				m->cmpR[curr_res] = m->g[gr].fstR + r;
				curr_res++;
			}
	return;
}

void output(BIGGEST_NUM *r, unsigned long nr, const char *filename) {
	FILE *f;
	unsigned long x;
	f = fopen(filename, "w");
	for (x = 0; x < nr; x++)
		fprintf(f, "%.18E\n", r[x]);
	fclose(f);
}

