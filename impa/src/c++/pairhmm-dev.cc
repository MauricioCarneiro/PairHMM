#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "timing.h"

#define MAX_TESTCASES_BUNCH_SIZE 100
#define VECTOR_SIZE 8

typedef struct
{
	int rslen, haplen, *q, *i, *d, *c;
	int *hap, *rs;
} testcase;

int normalize(char c)
{
	return (int)c - 33;
}

int read_testcase(testcase *tc)
{
	std::string hap, rs, q, i, d, c, i1, i2;

	if (!(std::cin >> hap >> rs >> q >> i >> d >> c >> i1 >> i2).good())
		return -1;

	tc->haplen = hap.size();
	tc->rslen = rs.size();

	int h = (tc->rslen+1);
	h = ((h + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;

	tc->hap = new int[tc->haplen + 2 * (h-1) + 1]();
	tc->hap += (h-1);
	tc->rs = new int[h+1]();
	tc->q = new int[h+1]();
	tc->i = new int[h+1]();
	tc->d = new int[h+1]();
	tc->c = new int[h+1]();

	for (int x = 0; x < tc->haplen; x++)
		tc->hap[tc->haplen-x-1] = hap[x];

	tc->hap += tc->haplen;

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
double compute_full_prob(testcase *tc, char *done);


template<>
double compute_full_prob<float>(testcase *tc, char *done)
{
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	/* constants */
	int sz = ((ROWS + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;

	float ph2pr[128], MM[sz + 1], GM[sz + 1], MX[sz + 1], XX[sz + 1], 
		MY[sz + 1], YY[sz + 1], pq[sz+1];
	for (int x = 0; x < 128; x++)
		ph2pr[x] = pow((float(10.0)), -(float(x)) / (float(10.0)));
	//	cell 0 of MM, GM, ... , YY is never used, since first row is just 
	//	"hard-coded" in the calculus (i.e.: not computed, just initialized).
	for (int r = 1; r < ROWS; r++)
	{
		int _i = (tc->i)[r-1] & 127;
		int _d = (tc->d)[r-1] & 127;
		int _c = (tc->c)[r-1] & 127;
		int _q = (tc->q)[r-1] & 127;
		//MM[r] = (float(1.0)) - ph2pr[(_i + _d) & 127];
		MM[r] = (float(1.0)) - ph2pr[_i]*ph2pr[_d];
		GM[r] = (float(1.0)) - ph2pr[_c];
		MX[r] = ph2pr[_i];
		XX[r] = ph2pr[_c];
		MY[r] = (r == ROWS - 1) ? (float(1.0)) : ph2pr[_d];
		YY[r] = (r == ROWS - 1) ? (float(1.0)) : ph2pr[_c];
		pq[r] = ph2pr[_q];
	}

	float M1[sz], M2[sz], M3[sz], *M, *Mp, *Mpp;
	float X1[sz], X2[sz], X3[sz], *X, *Xp, *Xpp;
	float Y1[sz], Y2[sz], Y3[sz], *Y, *Yp, *Ypp;
	Mpp = M1; Xpp = X1; Ypp = Y1;
	Mp = M2;  Xp = X2;  Yp = Y2;
	M = M3;   X = X3;   Y = Y3;


	/* first and second diagonals */
	float k = INITIAL_CONSTANT<float>() / (tc->haplen);

	Mpp[0] = (float(0.0));
	Xpp[0] = (float(0.0));
	Ypp[0] = k;
	Mp[0] = (float(0.0));
	Xp[0] = (float(0.0));
	Yp[0] = k;
	for (int r = 1; r < ROWS; r++)
	{
		Mpp[r] = (float(0.0));
		Xpp[r] = (float(0.0));
		Ypp[r] = (float(0.0));
		Mp[r] = (float(0.0));
		Xp[r] = (float(0.0));
		Yp[r] = (float(0.0));
	}

	/* main loop */
	float result = (float(0.0));
	for (int diag = 2; diag < (ROWS - 1) + COLS; diag++)
	{
		M[0] = (float(0.0));
		X[0] = (float(0.0));
		Y[0] = k;

		for (int rb = 1; rb < ROWS; rb += VECTOR_SIZE)
		{
			for (int r = rb; r < rb + VECTOR_SIZE; r++)
			{
				int _rs = tc->rs[r-1];
				int _hap = tc->hap[-diag+r];
				float distm = pq[r];
				if (_rs == _hap || _rs == 'N' || _hap == 'N')
					distm = (float(1.0)) - distm;

				M[r] = distm * (Mpp[r-1] * MM[r] + Xpp[r-1] * GM[r] + Ypp[r-1] * GM[r]);
				X[r] = Mp[r-1] * MX[r] + Xp[r-1] * XX[r];
				Y[r] = Mp[r] * MY[r] + Yp[r] * YY[r];
			}
		}

		result += M[ROWS-1] + X[ROWS-1];
		float *aux;
		aux = Mpp; Mpp = Mp; Mp = M; M = aux;
		aux = Xpp; Xpp = Xp; Xp = X; X = aux;
		aux = Ypp; Ypp = Yp; Yp = Y; Y = aux;
	}

	*done = (result > MIN_ACCEPTED<float>()) ? 1 : 0;

	return (double) (log10(result) - log10(INITIAL_CONSTANT<float>()));
}


template<>
double compute_full_prob<double>(testcase *tc, char *done)
{
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	/* constants */
	int sz = ((ROWS + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;

	double ph2pr[128], MM[sz + 1], GM[sz + 1], MX[sz + 1], XX[sz + 1], 
		MY[sz + 1], YY[sz + 1], pq[sz+1];
	for (int x = 0; x < 128; x++)
		ph2pr[x] = pow((double(10.0)), -(double(x)) / (double(10.0)));
	//	cell 0 of MM, GM, ... , YY is never used, since first row is just 
	//	"hard-coded" in the calculus (i.e.: not computed, just initialized).
	for (int r = 1; r < ROWS; r++)
	{
		int _i = (tc->i)[r-1] & 127;
		int _d = (tc->d)[r-1] & 127;
		int _c = (tc->c)[r-1] & 127;
		int _q = (tc->q)[r-1] & 127;
		//MM[r] = (double(1.0)) - ph2pr[(_i + _d) & 127];
		MM[r] = (double(1.0)) - ph2pr[_i]*ph2pr[_d];
		GM[r] = (double(1.0)) - ph2pr[_c];
		MX[r] = ph2pr[_i];
		XX[r] = ph2pr[_c];
		MY[r] = (r == ROWS - 1) ? (double(1.0)) : ph2pr[_d];
		YY[r] = (r == ROWS - 1) ? (double(1.0)) : ph2pr[_c];
		pq[r] = ph2pr[_q];
	}

	double M1[sz], M2[sz], M3[sz], *M, *Mp, *Mpp;
	double X1[sz], X2[sz], X3[sz], *X, *Xp, *Xpp;
	double Y1[sz], Y2[sz], Y3[sz], *Y, *Yp, *Ypp;
	Mpp = M1; Xpp = X1; Ypp = Y1;
	Mp = M2;  Xp = X2;  Yp = Y2;
	M = M3;   X = X3;   Y = Y3;


	/* first and second diagonals */
	double k = INITIAL_CONSTANT<double>() / (tc->haplen);

	Mpp[0] = (double(0.0));
	Xpp[0] = (double(0.0));
	Ypp[0] = k;
	Mp[0] = (double(0.0));
	Xp[0] = (double(0.0));
	Yp[0] = k;
	for (int r = 1; r < ROWS; r++)
	{
		Mpp[r] = (double(0.0));
		Xpp[r] = (double(0.0));
		Ypp[r] = (double(0.0));
		Mp[r] = (double(0.0));
		Xp[r] = (double(0.0));
		Yp[r] = (double(0.0));
	}

	/* main loop */
	double result = (double(0.0));
	for (int diag = 2; diag < (ROWS - 1) + COLS; diag++)
	{
		M[0] = (double(0.0));
		X[0] = (double(0.0));
		Y[0] = k;

		for (int rb = 1; rb < ROWS; rb += VECTOR_SIZE)
		{
			for (int r = rb; r < rb + VECTOR_SIZE; r++)
			{
				int _rs = tc->rs[r-1];
				int _hap = tc->hap[-diag+r];
				double distm = pq[r];
				if (_rs == _hap || _rs == 'N' || _hap == 'N')
					distm = (double(1.0)) - distm;

				M[r] = distm * (Mpp[r-1] * MM[r] + Xpp[r-1] * GM[r] + Ypp[r-1] * GM[r]);
				X[r] = Mp[r-1] * MX[r] + Xp[r-1] * XX[r];
				Y[r] = Mp[r] * MY[r] + Yp[r] * YY[r];
			}
		}

		result += M[ROWS-1] + X[ROWS-1];
		double *aux;
		aux = Mpp; Mpp = Mp; Mp = M; M = aux;
		aux = Xpp; Xpp = Xp; Xp = X; X = aux;
		aux = Ypp; Ypp = Yp; Yp = Y; Y = aux;
	}

	*done = (result > MIN_ACCEPTED<double>()) ? 1 : 0;

	return (double) (log10(result) - log10(INITIAL_CONSTANT<double>()));
}

int main()
{
	testcase tc[MAX_TESTCASES_BUNCH_SIZE];
	double result[MAX_TESTCASES_BUNCH_SIZE];
	char done[MAX_TESTCASES_BUNCH_SIZE];
	int num_tests;
    double kernel_time = 0.;
    Timing timer;

	do
	{
		num_tests = read_a_bunch_of_testcases(tc, MAX_TESTCASES_BUNCH_SIZE);
        timer.mark();
		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < num_tests; j++)
		{
			result[j] = compute_full_prob<float>(tc + j, done + j);
			if (!done[j])
				result[j] = compute_full_prob<double>(tc + j, done + j);
		}
        kernel_time += timer.elapsed();
		for (int j = 0; j < num_tests; j++)
			printf("%f\n", result[j]);
	} while (num_tests == MAX_TESTCASES_BUNCH_SIZE);
    fprintf(stderr, "%g\n", kernel_time);

	return 0;
}

