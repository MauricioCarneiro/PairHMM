#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "common.h"

void go1(Memory *mem, unsigned long i)
{
	unsigned long row, col, diag, R, H;
    double *p, *m, *mp, *mpp, *x, *xp, *xpp, *y, *yp, *ypp, *ph2pr, *swapper;
    char *ins, *del, *cont, *qual, *r, *h;

    H = mem->h[mem->cmpH[i]].H;
    h = mem->chunk + mem->h[mem->cmpH[i]].h;

    R = mem->r[mem->cmpR[i]].R;
    r = mem->chunk + mem->r[mem->cmpR[i]].r;
    qual = mem->chunk + mem->r[mem->cmpR[i]].qual;
    ins = mem->chunk + mem->r[mem->cmpR[i]].ins;
    del = mem->chunk + mem->r[mem->cmpR[i]].del;
    cont = mem->chunk + mem->r[mem->cmpR[i]].cont;

    ph2pr = mem->phred_to_prob;

    p = (double *) malloc((R+1)*15 * sizeof(double));
    m = p + (R+1)*6;
    mp = p + (R+1)*7;
    mpp = p + (R+1)*8;
    x = p + (R+1)*9;
    xp = p + (R+1)*10;
    xpp = p + (R+1)*11;
    y = p + (R+1)*12;
    yp = p + (R+1)*13;
    ypp = p + (R+1)*14;

    p[MtoM] = 1.0 - ph2pr[90];
    p[GapToM] = 1.0 - ph2pr[10];
    p[MtoX] = ph2pr[45];
    p[XtoX] = ph2pr[10];
    p[MtoY] = 1.0;
    p[YtoY] = 1.0;

    for (row = 1; row < R + 1; row++)
    {
        p[6*row+MtoM] = 1.0 - ph2pr[((ins[row-1] & 127) + (del[row-1] & 127)) & 127];
        p[6*row+GapToM] = 1.0 - ph2pr[cont[row-1] & 127];
        p[6*row+MtoX] = ph2pr[ins[row-1] & 127];
        p[6*row+XtoX] = ph2pr[cont[row-1] & 127];
        p[6*row+MtoY] = (row == R) ? 1.0 : ph2pr[del[row-1] & 127];
        p[6*row+YtoY] = (row == R) ? 1.0 : ph2pr[cont[row-1] & 127];
    }

    /*m[0] = 1.0;*/
	m[0] = pow(10, 300) / (H > R ? H - R + 1 : 1); 
    x[0] = 0.0;
    y[0] = 0.0;

    for (diag = 1; diag < R+H+1; diag++)
    {
        swapper = mpp; mpp = mp; mp = m; m = swapper;
        swapper = xpp; xpp = xp; xp = x; x = swapper;
        swapper = ypp; ypp = yp; yp = y; y = swapper;

        for (row = 0; row <= R; row++)
        {
            col = diag - row;
            if ((diag <= H + row) && (diag >= row))
            {
                m[row] = ((row == 0) || (col == 0)) ? 0.0 : \
                    ( \
						(r[row-1] == h[col-1] || r[row-1] == 'N' || h[col-1] == 'N') ? \
                        	1.0 - ph2pr[qual[row-1] & 127] : \
                        	ph2pr[qual[row-1] & 127] \
					) * (\
                        mpp[row-1] * p[6*row+MtoM] + \
                        xpp[row-1] * p[6*row+GapToM] + \
                        ypp[row-1] * p[6*row+GapToM] \
						);
                x[row] = (row == 0) ? 0.0 : mp[row-1] * p[6*row+MtoX] + xp[row-1] * p[6*row+XtoX];
                y[row] = (col == 0) ? 0.0 : mp[row] * p[6*row+MtoY] + yp[row] * p[6*row+YtoY];
            }
        }
    }

    /*mem->res[i] = log10(m[R] + x[R] + y[R]);*/
    mem->res[i] = log10(m[R] + x[R] + y[R]) - 300.0;

    free(p);

	return;
}

int main(int argc, char **argv)
{
    Memory m;
    unsigned long x;
	double t1, t2, t3, t4;
	
    if (argc != 3)
    {
        printf("\nUsage: <binary> <input file> <output cpu>\n\n");
        exit(0);
    }
	t1 = right_now();
    init_memory(argv[1], &m);
	t2 = right_now();
/* #pragma omp parallel for schedule(dynamic) */
	for (x = 0; x < m.nres; x++)
		go1(&m, x);
	t3 = right_now();
	output(m.res, m.nres, argv[2]);
	t4 = right_now();

	printf("INIT_MEMORY_TIME %e\n", 1000.0 * (t2 - t1));
	printf("COMPUTATION_TIME %e\n", 1000.0 * (t3 - t2));
	printf("OUTPUT_TIME %e\n", 1000.0 * (t4 - t3));
	printf("TOTAL_TIME %e\n", 1000.0 * (t4 - t1));
    return 0;
}

