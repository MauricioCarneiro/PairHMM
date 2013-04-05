#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <omp.h>

using std::cout;
using std::min;
using std::isprint;

#include "common.h"

typedef unsigned long ul;
typedef unsigned char uc;
typedef double dbl;

#define MAX_H 10240
#define MAX_CONCURRENT_COMPS 150

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
	ul buffstart, buffsz;
};

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void load_previous_results(ul j, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	ul tid = 0;//threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = 1;//blockDim.x * blockDim.y;
	ul c;

	for (c = ((j==0)?1:0) + tid; c < min(pv.H+1-j*DIAGS, DIAGS)+1; c+=nth)
	{
		pv.lastM[c] = pv.g_lastM[j*DIAGS+c-1];
		pv.lastX[c] = pv.g_lastX[j*DIAGS+c-1];
		pv.lastY[c] = pv.g_lastY[j*DIAGS+c-1];
	}
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void inline load_hap_data(ul j, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	ul tid = 0;//threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = 1;//blockDim.x * blockDim.y;
	ul c;

	for (c = tid; c < ROWS-1; c+=nth)
		pv.h[c] = pv.h[c + DIAGS];

	for (c=((j==0)?1:0)+tid; c < min(DIAGS, pv.H + 1 - j * DIAGS); c+=nth)
		pv.h[ROWS-1+c] = pv.chunk[pv.hap.h + j*DIAGS + c - 1];
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void inline load_read_data(ul i, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	ul tid = 0;//threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = 1;//blockDim.x * blockDim.y;
	ul r;

	for (r = ((i==0)?1:0) + tid; r < min(ROWS, pv.R+1-i*ROWS); r+=nth)
	{
		pv.r[r] = pv.chunk[pv.rs.r + i*ROWS + r - 1];
		pv.q[r] = pv.chunk[pv.rs.qual + i*ROWS + r - 1];
		pv.i[r] = pv.chunk[pv.rs.ins + i*ROWS + r - 1];
		pv.d[r] = pv.chunk[pv.rs.del + i*ROWS + r - 1];
		pv.c[r] = pv.chunk[pv.rs.cont + i*ROWS + r - 1];
	}
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void inline block(ul i, ul j, PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	ul tid = 0;//threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = 1;//blockDim.x * blockDim.y;
	dbl *sw, *m, *mp, *mpp, *x, *xp, *xpp, *y, *yp, *ypp;
	int c, r, nRows;
	ul diag;
	dbl mu, ml, mul, xu, xul, yl, yul, dist;
	dbl MX, XX, GapM, MM, MY, YY, Q;	
	ul k;

	(void) nth;
	m=pv.m; mp=pv.mp; mpp=pv.mpp;
	x=pv.x; xp=pv.xp; xpp=pv.xpp;
	y=pv.y; yp=pv.yp; ypp=pv.ypp;

	if (i > 0)
		load_previous_results<ROWS, DIAGS, BUFFSIZE>(j, pv);

	load_hap_data<ROWS, DIAGS, BUFFSIZE>(j, pv);

	nRows = min((int)ROWS, (int)(pv.R+1-i*ROWS));
	for (diag=0; diag < DIAGS; diag++) {
		for (r=0; r < nRows; r++) {
			c = diag - r + j * DIAGS;
			if (c < 0 || c > (int)pv.H) 
				continue;
			if (r == 0) { // First row of rowblock.
				if (i == 0) {	// First row of rowblock. First row of computation.
					if (c == 0) { // First row of rowblock. First row of computation. First column.
						/*m[0] = 1.0;*/
						m[r] = 1.0e300 / (pv.H > pv.R ? pv.H - pv.R + 1 : 1); 
						x[r] = 0.0;
						y[r] = 0.0;
					}
					else { // First row of rowblock. First row of computation. Not first column.
						m[r] = 0.0;
						x[r] = 0.0;
						y[r] = mp[0] + yp[0];
					}
				}
				else { // First row of rowblock. Not first row of computation.
					if (c == 0) { // First row of rowblock. Not first row of computation. First column.
						m[r] = 0.0;

						mu = pv.lastM[1];
						xu = pv.lastX[1];
						MX = pv.ph2pr[pv.i[r]];
						XX = pv.ph2pr[pv.c[r]];
						x[r] = mu * MX + xu * XX;

						y[r] = 0.0;
					}
					else { // First row of rowblock. Not first row of computation. Not first column.
						MM = 1.0 - pv.ph2pr[(pv.i[r] + pv.d[r]) & 127];
						GapM = 1.0 - pv.ph2pr[pv.c[r]];
						Q = pv.ph2pr[pv.q[r]];
						mul = pv.lastM[diag];
						xul = pv.lastX[diag];
						yul = pv.lastY[diag];
						dist = (pv.r[r]==pv.h[ROWS-1+diag-r] || pv.r[r]=='N' || pv.h[ROWS-1+diag-r]=='N') ? 1.0-Q : Q;
						m[r] = dist * (mul * MM + xul * GapM + yul * GapM);

						MX = pv.ph2pr[pv.i[r]];
						XX = pv.ph2pr[pv.c[r]];
						mu = pv.lastM[diag+1];
						xu = pv.lastX[diag+1];
						x[r] = mu * MX + xu * XX;

						MY = ((unsigned)r+i*ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.d[r]];
						YY = ((unsigned)r+i*ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.c[r]];
						ml = mp[r];
						yl = yp[r];
						y[r] = ml * MY + yl * YY;
					}
				}
			}
			else { // Not first row of rowblock.
				if (c == 0) { // Not first row of rowblock. First column.
					m[r] = 0.0;

					mu = mp[r-1];
					xu = xp[r-1];
					MX = pv.ph2pr[pv.i[r]];
					XX = pv.ph2pr[pv.c[r]];
					x[r] = mu * MX + xu * XX;

					y[r] = 0.0;
				}
				else { // Not first row of rowblock. Not first column.
					MM = 1.0 - pv.ph2pr[(pv.i[r] + pv.d[r]) & 127];
					GapM = 1.0 - pv.ph2pr[pv.c[r]];
					Q = pv.ph2pr[pv.q[r]];
					mul = mpp[r-1];
					xul = xpp[r-1];
					yul = ypp[r-1];
					dist = (pv.r[r]==pv.h[ROWS-1+diag-r] || pv.r[r]=='N' || pv.h[ROWS-1+diag-r]=='N') ? 1.0-Q : Q;
					m[r] = dist * (mul * MM + xul * GapM + yul * GapM);

					MX = pv.ph2pr[pv.i[r]];
					XX = pv.ph2pr[pv.c[r]];
					mu = mp[r-1];
					xu = xp[r-1];
					x[r] = mu * MX + xu * XX;

					MY = ((unsigned)r+i*ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.d[r]];
					YY = ((unsigned)r+i*ROWS == pv.R) ? 1.0 : pv.ph2pr[pv.c[r]];
					ml = mp[r];
					yl = yp[r];
					y[r] = ml * MY + yl * YY;
				}
			}

			if ((r == ROWS-1)&& (i < pv.nblockrows - 1))
			{
				pv.buffM[pv.buffsz] = m[r];
				pv.buffX[pv.buffsz] = x[r];
				pv.buffY[pv.buffsz] = y[r];
				pv.buffsz++;
			}
		}
		if (pv.buffsz==BUFFSIZE)
		{
			for (k=0; k<pv.buffsz; k++)
			{
				pv.g_lastM[pv.buffstart + k] = pv.buffM[k];
				pv.g_lastX[pv.buffstart + k] = pv.buffX[k];
				pv.g_lastY[pv.buffstart + k] = pv.buffY[k];
			}
			pv.buffstart+=BUFFSIZE;
			pv.buffsz=0;
		}

		sw=mpp; mpp=mp; mp=m; m=sw;
		sw=xpp; xpp=xp; xp=x; x=sw;
		sw=ypp; ypp=yp; yp=y; y=sw;
		if (tid==0)
		{
			pv.mpp=mpp; pv.mp=mp; pv.m=m;
			pv.xpp=xpp; pv.xp=xp; pv.x=x;
			pv.ypp=ypp; pv.yp=yp; pv.y=y;
		}
	}

	return;
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void test(Memory mem, ul compIndex, dbl *prevlines)
{
	ul tid = 0;//threadIdx.x + blockDim.x * threadIdx.y;
	ul nth = 1;//blockDim.x * blockDim.y;
	ul i, j, k;
	int id;
	dbl *lastM, *lastX, *lastY;

	id = omp_get_thread_num();

	lastM=prevlines + (id * 3 * MAX_H);
	lastX=lastM + MAX_H;
	lastY=lastX + MAX_H;

	PubVars<ROWS, DIAGS, BUFFSIZE> pv;

	for (i = 0+tid; i < 128; i+=nth)
		pv.ph2pr[i] = pow(10.0, -((double) i) / 10.0);

	ul cmph = mem.cmpH[compIndex];
	ul cmpr = mem.cmpR[compIndex];
	pv.rs = mem.r[cmpr];
	pv.hap = mem.h[cmph];
	pv.chunk = mem.chunk;
	pv.buffsz = 0;
	pv.buffstart = 0;
	pv.H = pv.hap.H;
	pv.R = pv.rs.R;
	pv.nblockrows = ((pv.R+1)+(ROWS-1)) / ROWS;
	pv.nblockcols = (((ROWS)+(pv.H+1)-1) + (DIAGS-1)) / DIAGS;
	pv.g_lastM = lastM;
	pv.g_lastX = lastX;
	pv.g_lastY = lastY;
	pv.m = pv.m1; pv.mp = pv.m2; pv.mpp = pv.m3;
	pv.y = pv.y1; pv.yp = pv.y2; pv.ypp = pv.y3;
	pv.x = pv.x1; pv.xp = pv.x2; pv.xpp = pv.x3;

	for (i=0; i < pv.nblockrows; i++)
	{
		load_read_data(i, pv);
		pv.buffstart=0;
		pv.buffsz=0;
		for (j=0; j < pv.nblockcols; j++)
			block<ROWS, DIAGS, BUFFSIZE>(i, j, pv);

		if (pv.buffsz > 0)
		{
			for (k=0; k<pv.buffsz; k++)
			{
				pv.g_lastM[pv.buffstart + k] = pv.buffM[k];
				pv.g_lastX[pv.buffstart + k] = pv.buffX[k];
				pv.g_lastY[pv.buffstart + k] = pv.buffY[k];
			}
		}
	}
	int resr = (int)pv.R - (int)(pv.nblockrows-1) * (int)ROWS;
	int swaps = ((int)DIAGS * (int)pv.nblockcols -resr - (int)(pv.H + 1)) % 3;
	if (swaps == 2)
		mem.res[compIndex] = log10(pv.m[resr] + pv.x[resr] + pv.y[resr]) - 300.0;
	else if (swaps == 0)
		mem.res[compIndex] = log10(pv.mp[resr] + pv.xp[resr] + pv.yp[resr]) - 300.0;
	else
		mem.res[compIndex] = log10(pv.mpp[resr] + pv.xpp[resr] + pv.ypp[resr]) - 300.0;
	return;
}

int main(int argc, char **argv)
{
    Memory m;
	ul x;
	dbl *lastlines = (dbl *) malloc(MAX_CONCURRENT_COMPS * 3 * sizeof(dbl) * MAX_H);

    if (argc != 3)
    {
        printf("\nUsage: <binary> <input file> <output>\n\n");
        exit(0);
    }
	init_memory(argv[1], &m);

#pragma omp parallel for schedule(dynamic) 
	for (x = 0; x < m.nres; x++)
		test<128, 300, 30>(m, x, lastlines);

	output(m.res, m.nres, argv[2]);

	free(lastlines);
    return 0;
}

/*
template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void debug_read_data(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	cout << "r = [";
	for (ul x = 0; x < ROWS; x++)
		cout << pv.r[x];
	cout << "]\n";
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void debug_buffers_data(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	for (ul x = 0; x < BUFFSIZE; x++) printf("%e ", pv.buffM[x]); printf("\n");
	for (ul x = 0; x < BUFFSIZE; x++) printf("%e ", pv.buffX[x]); printf("\n");
	for (ul x = 0; x < BUFFSIZE; x++) printf("%e ", pv.buffY[x]); printf("\n");
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void debug_prev_data(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	for (ul x = 0; x < pv.H+1; x++) printf("%e ", pv.g_lastM[x]); printf("\n");
	for (ul x = 0; x < pv.H+1; x++) printf("%e ", pv.g_lastX[x]); printf("\n");
	for (ul x = 0; x < pv.H+1; x++) printf("%e ", pv.g_lastY[x]); printf("\n");
}

template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void debug_prev_data_sh(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	for (ul x = 0; x < DIAGS+1; x++) printf("%e ", pv.lastM[x]); printf("\n");
	for (ul x = 0; x < DIAGS+1; x++) printf("%e ", pv.lastX[x]); printf("\n");
	for (ul x = 0; x < DIAGS+1; x++) printf("%e ", pv.lastY[x]); printf("\n");
}


template <ul ROWS, ul DIAGS, ul BUFFSIZE>
void debug_hap_data(PubVars<ROWS, DIAGS, BUFFSIZE> &pv)
{
	cout << "hap = [";
	for (ul x = 0; x < ROWS+DIAGS-1; x++)
		cout << pv.h[x];
	cout << "]\n";
}
*/

/*
#define DBG_ROWS 39
#define DBG_COLS 246

double DEBUGm[DBG_ROWS][DBG_COLS];
double DEBUGx[DBG_ROWS][DBG_COLS];
double DEBUGy[DBG_ROWS][DBG_COLS];

void DEBUG_PRINT_MATRICES()
{
	FILE *fm, *fx, *fy;
	int r, c;

	fm = fopen("M.txt", "w");
	fx = fopen("X.txt", "w");
	fy = fopen("Y.txt", "w");

	for (r = 0; r < DBG_ROWS; r++)
	{
		for (c = 0; c < DBG_COLS; c++)
		{
			fprintf(fm, "%e ", DEBUGm[r][c]);
			fprintf(fx, "%e ", DEBUGx[r][c]);
			fprintf(fy, "%e ", DEBUGy[r][c]);
		}
		fprintf(fm, "\n");
		fprintf(fx, "\n");
		fprintf(fy, "\n");
	}

	fclose(fm);
	fclose(fx);
	fclose(fy);
}

void SET_VALUES_DEBUG_MATRICES(int r, int c, double m, double x, double y)
{
	DEBUGm[r][c] = m;
	DEBUGx[r][c] = x;
	DEBUGy[r][c] = y;
}
*/

