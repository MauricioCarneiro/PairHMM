#ifndef COMMON_H
#define COMMON_H

/*
#define REAL_NUMBER double
#define CONST(v) ((REAL_NUMBER) v)
#define KERNEL_POW(u,v) pow(u,v)
#define KERNEL_LOG10(v) log10(v)
#define INITIAL_CONSTANT (ldexp(1.0, 1020.0))
#define LOG10_INITIAL_CONSTANT (log10(INITIAL_CONSTANT))
*/

#define REAL_NUMBER float
#define CONST(v) ((REAL_NUMBER) v)
#define KERNEL_POW(u,v) powf(u,v)
#define KERNEL_LOG10(v) log10f(v)
#define INITIAL_CONSTANT (ldexpf(1.0, 120))
#define LOG10_INITIAL_CONSTANT (log10f(INITIAL_CONSTANT))

typedef struct
{
	unsigned long nR, nH, fstR, fstH;
} Group;

typedef struct
{
	unsigned long R;
	int r, qual, ins, del, cont;
} ReadSequence;

typedef struct
{
	unsigned long H;
	int h;
} Haplotype;

typedef struct
{
	Group *g;
	Haplotype *h;
	ReadSequence *r;
	char *chunk;
    unsigned long *cmpH;
    unsigned long *cmpR;
    REAL_NUMBER *res;
    REAL_NUMBER *phred_to_prob;
    unsigned long ng, nh, nr, chunk_sz, nres;
} Memory;

double right_now(void);
int ndigs(char c);
int read_numbers(const char *str, int n, char *dest);
void init_memory(const char *fn, Memory *m);
void show_time(const char *txt);
void output(REAL_NUMBER *r, unsigned long nr, const char *filename);

#endif
