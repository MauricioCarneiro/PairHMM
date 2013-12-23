#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>


#include "input.h"

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
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	/* constants */
	int sz = ((ROWS + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;

	NUMBER ph2pr[128], MM[sz + 1], GM[sz + 1], MX[sz + 1], XX[sz + 1], 
		MY[sz + 1], YY[sz + 1], pq[sz+1];
	for (int x = 0; x < 128; x++)
		ph2pr[x] = pow((NUMBER(10.0)), -(NUMBER(x)) / (NUMBER(10.0)));
	//	cell 0 of MM, GM, ... , YY is never used, since first row is just 
	//	"hard-coded" in the calculus (i.e.: not computed, just initialized).
	for (int r = 1; r < ROWS; r++)
	{
		int _i = (tc->i)[r-1] & 127;
		int _d = (tc->d)[r-1] & 127;
		int _c = (tc->c)[r-1] & 127;
		int _q = (tc->q)[r-1] & 127;
		//MM[r] = (NUMBER(1.0)) - ph2pr[(_i + _d) & 127];
		MM[r] = (NUMBER(1.0)) - ph2pr[_i]*ph2pr[_d];
		GM[r] = (NUMBER(1.0)) - ph2pr[_c];
		MX[r] = ph2pr[_i];
		XX[r] = ph2pr[_c];
		MY[r] = (r == ROWS - 1) ? (NUMBER(1.0)) : ph2pr[_d];
		YY[r] = (r == ROWS - 1) ? (NUMBER(1.0)) : ph2pr[_c];
		pq[r] = ph2pr[_q];
	}

	NUMBER M1[sz], M2[sz], M3[sz], *M, *Mp, *Mpp;
	NUMBER X1[sz], X2[sz], X3[sz], *X, *Xp, *Xpp;
	NUMBER Y1[sz], Y2[sz], Y3[sz], *Y, *Yp, *Ypp;
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
	for (int r = 1; r < ROWS; r++)
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
	for (int diag = 2; diag < (ROWS - 1) + COLS; diag++)
	{
		M[0] = (NUMBER(0.0));
		X[0] = (NUMBER(0.0));
		Y[0] = k;

		for (int rb = 1; rb < ROWS; rb += VECTOR_SIZE)
		{
			for (int r = rb; r < rb + VECTOR_SIZE; r++)
			{
				int _rs = tc->rs[r-1];
				int _hap = tc->hap[diag-r-1];
				NUMBER distm = pq[r];
				if (_rs == _hap || _rs == 'N' || _hap == 'N')
					distm = (NUMBER(1.0)) - distm;

				M[r] = distm * (Mpp[r-1] * MM[r] + Xpp[r-1] * GM[r] + Ypp[r-1] * GM[r]);
				X[r] = Mp[r-1] * MX[r] + Xp[r-1] * XX[r];
				Y[r] = Mp[r] * MY[r] + Yp[r] * YY[r];
			}
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

