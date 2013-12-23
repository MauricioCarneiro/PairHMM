#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "input.h"

#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5

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

	NUMBER ph2pr[128];
	for (int x = 0; x < 128; x++)
		ph2pr[x] = pow((NUMBER(10.0)), -(NUMBER(x)) / (NUMBER(10.0)));

	NUMBER M[ROWS][COLS + 2 * ROWS - 2];
	NUMBER X[ROWS][COLS + 2 * ROWS - 2];
	NUMBER Y[ROWS][COLS + 2 * ROWS - 2];
	NUMBER p[ROWS][6];

	p[0][MM] = (NUMBER(0.0));
	p[0][GapM] = (NUMBER(0.0));
	p[0][MX] = (NUMBER(0.0));
	p[0][XX] = (NUMBER(0.0));
	p[0][MY] = (NUMBER(0.0));
	p[0][YY] = (NUMBER(0.0));
	for (r = 1; r < ROWS; r++)
	{
		int _i = tc->i[r-1] & 127;
		int _d = tc->d[r-1] & 127;
		int _c = tc->c[r-1] & 127;
		p[r][MM] = (NUMBER(1.0)) - ph2pr[(_i + _d) & 127];
		p[r][GapM] = (NUMBER(1.0)) - ph2pr[_c];
		p[r][MX] = ph2pr[_i];
		p[r][XX] = ph2pr[_c];
		p[r][MY] = (r == ROWS - 1) ? (NUMBER(1.0)) : ph2pr[_d];
		p[r][YY] = (r == ROWS - 1) ? (NUMBER(1.0)) : ph2pr[_c];
	}

	/* first row */
	NUMBER k = INITIAL_CONSTANT<NUMBER>() / (tc->haplen);
	for (c = 0 + ROWS - 1; c < COLS + ROWS - 1; c++)
	{
		M[0][c] = (NUMBER(0.0));
		X[0][c] = (NUMBER(0.0));
		Y[0][c] = k;
	}

	/* first and second diagonal */
	for (r = 1; r < ROWS; r++)
	{
		M[r][(ROWS-1)-r] = (NUMBER(0.0));
		M[r][(ROWS-1)-r+1] = (NUMBER(0.0));
		X[r][(ROWS-1)-r] = (NUMBER(0.0));
		X[r][(ROWS-1)-r+1] = (NUMBER(0.0));
		Y[r][(ROWS-1)-r] = (NUMBER(0.0));
		Y[r][(ROWS-1)-r+1] = (NUMBER(0.0));
	}

	for (diag = 2; diag < (ROWS-1)+COLS; diag++)
		for (r = 1; r < ROWS; r++)
		{
			c = (ROWS-1)+diag-r;

			char _rs = tc->rs[r-1];
			char _hap = tc->hap[c-1-(ROWS-1)];
			int _q = tc->q[r-1] & 127;
			NUMBER distm = ph2pr[_q];
			if (_rs == _hap || _rs == 'N' || _hap == 'N')
				distm = (NUMBER(1.0)) - distm;
			M[r][c] = distm * (M[r-1][c-1] * p[r][MM] + X[r-1][c-1] * p[r][GapM] + Y[r-1][c-1] * p[r][GapM]);
			X[r][c] = M[r-1][c] * p[r][MX] + X[r-1][c] * p[r][XX];
			Y[r][c] = M[r][c-1] * p[r][MY] + Y[r][c-1] * p[r][YY];
		}

	NUMBER result = (NUMBER(0.0));
	for (c = 0 + (ROWS-1); c < COLS + (ROWS-1); c++)
		result += M[ROWS-1][c] + X[ROWS-1][c];

	/* DEBUG */
	#ifdef DEBUG_MODE
	printf(" M \n---\n");
	for (r = 0; r < ROWS; r++) {
		for (c = 0; c < COLS + 2 * ROWS - 1; c++)
			if ((ROWS-1-r > c) || (c >= COLS + 2 * ROWS - 2 - r))
				printf("          ");
			else
				printf(" %09.2e", M[r][c]);
		printf("\n");
	}
	printf(" X \n---\n");
	for (r = 0; r < ROWS; r++) {
		for (c = 0; c < COLS + 2 * ROWS - 1; c++)
			if ((ROWS-1-r > c) || (c >= COLS + 2 * ROWS - 2 - r))
				printf("          ");
			else
				printf(" %09.2e", X[r][c]);
		printf("\n");
	}
	printf(" Y \n---\n");
	for (r = 0; r < ROWS; r++) {
		for (c = 0; c < COLS + 2 * ROWS - 1; c++)
			if ((ROWS-1-r > c) || (c >= COLS + 2 * ROWS - 2 - r))
				printf("          ");
			else
				printf(" %09.2e", Y[r][c]);
		printf("\n");
	}
	#endif

	*done = (result > MIN_ACCEPTED<NUMBER>()) ? 1 : 0;

	return (double) (log10(result) - log10(INITIAL_CONSTANT<NUMBER>()));
}

int main()
{
	testcase tc[MAX_TESTCASES_BUNCH_SIZE];
	double result[MAX_TESTCASES_BUNCH_SIZE];
	char done[MAX_TESTCASES_BUNCH_SIZE];
	int num_tests;

	#ifdef DEBUG_MODE
		num_tests = read_a_bunch_of_testcases(tc, MAX_TESTCASES_BUNCH_SIZE);
		(void)num_tests;
		result[0] = compute_full_prob<float>(tc + 0, done + 0);
		if (!done[0])
			result[0] = compute_full_prob<double>(tc + 0, done + 0);
		printf("%f\n", result[0]);
	#else
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
	#endif

	return 0;
}

