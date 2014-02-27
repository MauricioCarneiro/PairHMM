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

//#define DEBUG_MODE
//#define DEBUG_ROWS

#include <vector>
using std::vector;


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

#ifdef DEBUG_ROWS
vector<int> rows_to_debug;
rows_to_debug.push_back(0);
rows_to_debug.push_back(1);
rows_to_debug.push_back(2);
rows_to_debug.push_back(32);
rows_to_debug.push_back(33);
rows_to_debug.push_back(34);
NUMBER debug[rows_to_debug.size()][10];
for (int y = 0; y < rows_to_debug.size(); y++) for (int x = 0; x < 10; x++) debug[y][x] = 5.0;
#endif
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

#ifdef DEBUG_ROWS
int _r = r, _c = c-(ROWS-1);
int thisindex = -1;
for (int x = 0; x < rows_to_debug.size(); x++)
	if (_r == rows_to_debug[x])
		thisindex = x;
if ((thisindex != -1) && (_c >= 0 && _c < 8))
	debug[thisindex][_c] = M[r][c];
#endif

#if 0 //ndef DEBUG_ROWS
int _r = r, _c = c-(ROWS-1);
if ((_r == 32 && _c == 1) || (_r == 33 && _c == 2))
{
printf("[%d][%d]\n", r, _c);
printf("m=%09.2e\n", M[r][c]);
printf("x=%09.2e\n", X[r][c]);
printf("y=%09.2e\n", Y[r][c]);
printf("rs=%c\n", _rs);
printf("q=%d. COEF=%09.2e\n", _q, distm);
printf("ph2pr[%d]=%09.2e\n", _q, ph2pr[_q]);
printf("hap=%c\n", _hap);
printf("mm=%09.2e\n", p[r][MM]);
printf("gm=%09.2e\n", p[r][GapM]);
printf("mx=%09.2e\n", p[r][MX]);
printf("xx=%09.2e\n", p[r][XX]);
printf("my=%09.2e\n", p[r][MY]);
printf("yy=%09.2e\n", p[r][YY]);
printf("M[%d][%d]=%09.2e\n", r-1, _c-1, M[r-1][c-1]);
printf("X[%d][%d]=%09.2e\n", r-1, _c-1, X[r-1][c-1]);
printf("Y[%d][%d]=%09.2e\n", r-1, _c-1, Y[r-1][c-1]);
printf("\n\n");
}
#endif
		}

#ifdef DEBUG_ROWS
for (int y = 0; y < rows_to_debug.size(); y++) {
	printf("[%2d]    ", rows_to_debug[y]);
	for (int x = 0; x < 8; x++)
		if (debug[y][x] != 5.0)
			printf("%09.2e    ", debug[y][x]);
		else
			printf("             ");
	printf("\n");
}
#endif

	NUMBER result = (NUMBER(0.0));
	for (c = 0 + (ROWS-1); c < COLS + (ROWS-1); c++)
		result += M[ROWS-1][c] + X[ROWS-1][c];

	#if 0
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

