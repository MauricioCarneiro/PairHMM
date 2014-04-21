#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "input.h"
#include <iostream>
#include <vector>

int normalize(char c)
{
	return ((int) (c - 33));
}

int readv2(testcase *tc, std::istream& infile, int n_read, int n_hap) 
/* reads from infile reads and haps in the following format
#reads   #haps
rs q i d c
rs q i d c
...
rs q i d c
h 
h
...
h
<eof>

To read a second set like above, call readv2 again (with new n_read, n_hap)
*/
{
    static int read_inx=0;
    static bool init = false;
    static int remaining;
    static std::vector<testcase*>cases;
    static std::string hap;
    std::string rs, q, i, d, c;
    if (!init) 
    {
       remaining = n_read * n_hap;
       for (int z=0;z<n_read;z++)
       {
           testcase* this_tc = new testcase;
           if (!(infile >> rs >> q >> i >> d >> c).good())
               return -1;
           this_tc->rslen = rs.size();

       	  int sz = 1 + ((this_tc->rslen + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;

           this_tc->rs = new char[sz]();
           this_tc->q = new int[sz]();
           this_tc->i = new int[sz]();
           this_tc->d = new int[sz]();
           this_tc->c = new int[sz]();


           for (int x = 0; x < this_tc->rslen; x++)
           {
               this_tc->rs[x] = rs[x];
               this_tc->q[x] = normalize((q.c_str())[x]);
               this_tc->q[x] = this_tc->q[x] < 6 ? 6 : (this_tc->q[x]) & 127;
               this_tc->i[x] = normalize((i.c_str())[x]);
               this_tc->d[x] = normalize((d.c_str())[x]);
               this_tc->c[x] = normalize((c.c_str())[x]);
           }
           this_tc->hap = 0;
           this_tc->haplen = 0;
           cases.push_back(this_tc);
       }
       init = true;
       read_inx = n_read; //to trigger new hap read
    }
    if (remaining==0) {
       for (int z=0; z<n_read;z++) delete(cases[z]);
       init = false;
       read_inx=0;
       cases.clear();
       return -1;
    }
    
    if (read_inx == n_read) 
    {
       if (!(infile >> hap).good())
       {
           for (int z=0; z<n_read;z++) delete(cases[z]);
           return -1;
       }
       read_inx=0;
    }
    //if (cases[read_inx]->hap) printf("cases[%d]->hap = %s\n", read_inx, cases[read_inx]->hap);
    //else printf("cases[%d]->hap = (null)\n", read_inx);
    *tc = *cases[read_inx++]; //copy constructor
    tc->haplen = hap.size();
	 int sz = 1 + ((tc->rslen + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;
    if (tc->hap) {
       delete [] tc->hap;
    }
    tc->hap = new char[tc->haplen + 2 * (sz-1) + 1]();
    tc->hap += (sz-1);
    for (int x = 0; x < tc->haplen; x++)
        tc->hap[x] = hap[x];
    
    remaining--;
    return 0;
}

int read_testcase(testcase *tc, std::istream& infile, bool in_new)
{
    static bool newfile = true;
    static int version=0;
    static int n_read, n_hap;
    
    std::string hap, rs, q, i, d, c, i1, i2;

    if (in_new) newfile = true;
    if (newfile) {
       newfile = false;
       if (!(infile >> hap >> rs).good())
          return -1;
       if (hap == "#Version" && rs=="2") {
          if (!(infile >> hap >> rs).good())
             return -1;
          n_read = atoi(hap.c_str()); 
          n_hap = atoi(rs.c_str()); 
          if (n_read <=0) return -1;
          version = 2;
          return readv2(tc, infile, n_read, n_hap);
       } else {
          if (!(infile >> q >> i >> d >> c).good())
             return -1;
       }
    } else if (version==2) {
       if (-1 == readv2(tc, infile, n_read, n_hap)) 
       {
          if (!(infile >> hap >> rs).good()) return -1;
          n_read = atoi(hap.c_str());
          n_hap = atoi(rs.c_str());
          return readv2(tc, infile, n_read, n_hap);
       } else {
          //printf("tc = %p\n", tc);
          //tc->display();
          return 0;
       }
    } else {
       if (!(infile >> hap >> rs >> q >> i >> d >> c).good())
           return -1;
    }

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

int read_a_bunch_of_testcases(testcase *tc, int max_bunch_size, std::istream& infile)
{
	int num_tests = 0;
	for (num_tests = 0; 
		(num_tests < max_bunch_size) && (read_testcase(tc + num_tests, infile) == 0); 
		num_tests++);
	return num_tests;
}

