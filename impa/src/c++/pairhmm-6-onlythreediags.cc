#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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
	char *q, *i, *d, *c, *line = NULL;
	int _q, _i, _d, _c;
	int x, size = 0;
	ssize_t read;

	read = getline(&line, (size_t *) &size, stdin);
	if (read == -1)
		return -1;

	tc->hap = (char *) malloc(size);
	tc->rs = (char *) malloc(size);
	q = (char *) malloc(size);
	i = (char *) malloc(size);
	d = (char *) malloc(size);
	c = (char *) malloc(size);

	if (sscanf(line, "%s %s %s %s %s %s\n", tc->hap, tc->rs, q, i, d, c) != 6)
		return -1;

	tc->haplen = strlen(tc->hap);
	tc->rslen = strlen(tc->rs);
	tc->q = (int *) malloc(sizeof(int) * tc->rslen);
	tc->i = (int *) malloc(sizeof(int) * tc->rslen);
	tc->d = (int *) malloc(sizeof(int) * tc->rslen);
	tc->c = (int *) malloc(sizeof(int) * tc->rslen);

	for (x = 0; x < tc->rslen; x++)
	{
		_q = normalize(q[x]);
		_i = normalize(i[x]);
		_d = normalize(d[x]);
		_c = normalize(c[x]);
		tc->q[x] = (_q < 6) ? 6 : _q;
		tc->i[x] = _i;
		tc->d[x] = _d;
		tc->c[x] = _c;
	}

	free(q);
	free(i);
	free(d);
	free(c);
	free(line);

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

