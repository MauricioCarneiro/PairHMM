#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "input.h"

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

	int sz = 1 + ((tc->rslen + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;

    tc->hap = new char[tc->haplen + 2 * (sz-1) + 1]();
    tc->hap += (sz-1);
    tc->rs = new char[sz]();
    tc->q = new int[sz]();
    tc->i = new int[sz]();
    tc->d = new int[sz]();
    tc->c = new int[sz]();

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

int read_a_bunch_of_testcases(testcase *tc, int max_bunch_size)
{
	int num_tests = 0;
	for (num_tests = 0; 
		(num_tests < max_bunch_size) && (read_testcase(tc + num_tests) == 0); 
		num_tests++);
	return num_tests;
}

