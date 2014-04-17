#ifndef INPUT_H
#define INPUT_H

#define MAX_TESTCASES_BUNCH_SIZE 100

#define VECTOR_SIZE 8

#include <fstream>

struct testcase
{
	int rslen, haplen, *q, *i, *d, *c;
   double foo, bar, baz, boop;
	char *hap, *rs;
   void display()
   {
      printf("hap (len %d) = %s\n", haplen, hap);
      printf("rs (len %d) = %s\n", rslen, rs);
      printf("q = ");
      for (int z=0;z<rslen;z++) printf("%d  ", q[z]);
      printf("\n");
      printf("i = ");
      for (int z=0;z<rslen;z++) printf("%d  ", i[z]);
      printf("\n");
      printf("d = ");
      for (int z=0;z<rslen;z++) printf("%d  ", d[z]);
      printf("\n");
      printf("c = ");
      for (int z=0;z<rslen;z++) printf("%d  ", c[z]);
      printf("\n");
   }
   testcase& operator=(const testcase& src)
   {
       rslen = src.rslen;
       haplen = src.haplen;
	    int sz = 1 + ((rslen + VECTOR_SIZE - 1) / VECTOR_SIZE) * VECTOR_SIZE;
       if (rslen > 0) 
       {
          q = new int[sz]();
          i = new int[sz]();
          c = new int[sz]();
          d = new int[sz]();
       }
       if (src.hap) 
       {
          hap = new char[haplen + 2 * (sz-1) + 1]();
          hap += (sz-1);
          strncpy(hap, src.hap, src.haplen);
       } else {
          hap = 0;
       }
       if (src.rs) 
       {
          rs = new char[sz]();
          strncpy(rs, src.rs, src.rslen);
       } else {
          rs = 0;
       }
       for (int z=0;z<rslen;z++) 
       {
           q[z] = src.q[z];
           i[z] = src.i[z];
           d[z] = src.d[z];
           c[z] = src.c[z];
       }
       return *this;
   } 
   testcase() {
       hap = rs = 0;
       q = i = d = c = 0;
       haplen = rslen = 0;
   }
#if 0
   ~testcase() {
       printf("deleting %p. ", this);
       if (hap) printf("Freeing %p\n", hap);
       else printf("hap = %p\n", hap);
       if (hap) delete(hap); 
       if (rs) delete(rs);
       if (q) delete(q);
       if (i) delete(i);
       if (d) delete(d);
       if (c) delete(c);
       hap = rs = 0;
       q = i = d = c = 0;
   }
#endif
};

int read_testcase(testcase *, std::istream&, bool=false);
int read_a_bunch_of_testcases(testcase *, int, std::fstream&);

#endif

