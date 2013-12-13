#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "timing.h"

/*
#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5
*/

#define MAX_TESTCASES_BUNCH_SIZE 100

/*
	q: read quality
	i: insertion penalty
	d: deletion penalty
	c: gap continuation penalty
*/

typedef struct 
{
	int rslen, haplen, *q, *i, *d, *c;
	char *hap, *rs;
} testcase;

int normalize(char c)
{
	return ((int) (c - 33));
}

int read_testcase(testcase *tc)
{
    std::string hap, rs, q, i, d, c, i1, i2;

    if (!(std::cin >> hap >> rs >> q >> i >> d >> c >> i1 >> i2).good())
        return -1;

    tc->haplen = hap.size();
    tc->rslen = rs.size();

    int h = (tc->rslen+1);

    tc->hap = new char[tc->haplen + 2 * (h-1) + 1]();
    tc->hap += (h-1);
    tc->rs = new char[h+1]();
    tc->q = new int[h+1]();
    tc->i = new int[h+1]();
    tc->d = new int[h+1]();
    tc->c = new int[h+1]();

    for (int x = 0; x < tc->haplen; x++)
        tc->hap[x] = hap[x];

    for (int x = 0; x < tc->rslen; x++)
    {
        tc->rs[x] = rs[x];
        tc->q[x] = normalize((q.c_str())[x]);
        tc->q[x] = tc->q[x] < 6 ? 6 : (tc->q[x]) & 127;
        tc->i[x] = normalize((i.c_str())[x]);
        tc->d[x] = normalize((d.c_str())[x]);
        tc->c[x] = normalize((c.c_str())[x]);
    }

	return 0;
}

inline int read_a_bunch_of_testcases(testcase *tc, int max_bunch_size)
{
	int num_tests = 0;
	for (num_tests = 0; 
		(num_tests < max_bunch_size) && (read_testcase(tc + num_tests) == 0); 
		num_tests++);
	return num_tests;
}

template<class T>
inline T INITIAL_CONSTANT();

template<>
inline float INITIAL_CONSTANT<float>()
{
	return 1e32f;
}

template<>
inline double INITIAL_CONSTANT<double>()
{
	return ldexp(1.0, 1020);
}

template<class T>
inline T MIN_ACCEPTED();

template<>
inline float MIN_ACCEPTED<float>()
{
	return 1e-28f;
}

template<>
inline double MIN_ACCEPTED<double>()
{
	return 0.0;
}

template<class NUMBER>
double compute_full_prob(testcase *tc, char *done)
{
	int r, c, diag;
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	/* constants */
	NUMBER ph2pr[128], MM[ROWS], GM[ROWS], MX[ROWS], XX[ROWS], MY[ROWS], YY[ROWS];
	for (int x = 0; x < 128; x++)
		ph2pr[x] = pow((NUMBER(10.0)), -(NUMBER(x)) / (NUMBER(10.0)));
	//	cell 0 of MM, GM, ... , YY is never used, since first row is just 
	//	"hard-coded" in the calculus (i.e.: not computed, just initialized).
	for (r = 1; r < ROWS; r++)
	{
		int _i = tc->i[r-1] & 127;
		int _d = tc->d[r-1] & 127;
		int _c = tc->c[r-1] & 127;
		MM[r] = (NUMBER(1.0)) - ph2pr[(_i + _d) & 127];
		GM[r] = (NUMBER(1.0)) - ph2pr[_c];
		MX[r] = ph2pr[_i];
		XX[r] = ph2pr[_c];
		MY[r] = (r == ROWS - 1) ? (NUMBER(1.0)) : ph2pr[_d];
		YY[r] = (r == ROWS - 1) ? (NUMBER(1.0)) : ph2pr[_c];
	}

	NUMBER M1[ROWS], M2[ROWS], M3[ROWS], *M, *Mp, *Mpp;
	NUMBER X1[ROWS], X2[ROWS], X3[ROWS], *X, *Xp, *Xpp;
	NUMBER Y1[ROWS], Y2[ROWS], Y3[ROWS], *Y, *Yp, *Ypp;
	Mpp = M1; Xpp = X1; Ypp = Y1;
	Mp = M2;  Xp = X2;  Yp = Y2;
	M = M3;   X = X3;   Y = Y3;

	/* first and second diagonals */
	NUMBER k = INITIAL_CONSTANT<NUMBER>() / (tc->haplen);

	Mpp[0] = (NUMBER(0.0));
	Xpp[0] = (NUMBER(0.0));
	Ypp[0] = k;
	Mp[0] = (NUMBER(0.0));
	Xp[0] = (NUMBER(0.0));
	Yp[0] = k;
	for (r = 1; r < ROWS; r++)
	{
		Mpp[r] = (NUMBER(0.0));
		Xpp[r] = (NUMBER(0.0));
		Ypp[r] = (NUMBER(0.0));
		Mp[r] = (NUMBER(0.0));
		Xp[r] = (NUMBER(0.0));
		Yp[r] = (NUMBER(0.0));
	}

	/* main loop */
	NUMBER result = (NUMBER(0.0));
	for (diag = 2; diag < (ROWS-1)+COLS; diag++)
	{
		M[0] = (NUMBER(0.0));
		X[0] = (NUMBER(0.0));
		Y[0] = k;

		for (r = 1; r < ROWS; r++)
		{
			c = (ROWS-1)+diag-r;

			char _rs = tc->rs[r-1];
			char _hap = tc->hap[c-1-(ROWS-1)];
			int _q = tc->q[r-1] & 127;
			NUMBER distm = ph2pr[_q];
			if (_rs == _hap || _rs == 'N' || _hap == 'N')
				distm = (NUMBER(1.0)) - distm;

			M[r] = distm * (Mpp[r-1] * MM[r] + Xpp[r-1] * GM[r] + Ypp[r-1] * GM[r]);
			X[r] = Mp[r-1] * MX[r] + Xp[r-1] * XX[r];
			Y[r] = Mp[r] * MY[r] + Yp[r] * YY[r];
		}

		result += M[ROWS-1] + X[ROWS-1];
		NUMBER *aux;
		aux = Mpp; Mpp = Mp; Mp = M; M = aux;
		aux = Xpp; Xpp = Xp; Xp = X; X = aux;
		aux = Ypp; Ypp = Yp; Yp = Y; Y = aux;
	}

	*done = (result > MIN_ACCEPTED<NUMBER>()) ? 1 : 0;

	return (double) (log10(result) - log10(INITIAL_CONSTANT<NUMBER>()));
}

int main()
{
    Timing t;
	testcase tc[MAX_TESTCASES_BUNCH_SIZE];
	double result[MAX_TESTCASES_BUNCH_SIZE];
	char done[MAX_TESTCASES_BUNCH_SIZE];
	int num_tests;

	do
	{
		num_tests = read_a_bunch_of_testcases(tc, MAX_TESTCASES_BUNCH_SIZE);

		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < num_tests; j++)
		{
			result[j] = compute_full_prob<float>(tc + j, done + j);
			if (!done[j])
				result[j] = compute_full_prob<double>(tc + j, done + j);
		}
		for (int j = 0; j < num_tests; j++)
			printf("%f\n", result[j]);
	} while (num_tests == MAX_TESTCASES_BUNCH_SIZE);

	return 0;
}

