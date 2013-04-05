#ifndef COMMON_H
#define COMMON_H

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
    double *res;
    double *phred_to_prob;
    unsigned long ng, nh, nr, chunk_sz, nres;
} Memory;

double right_now(void);
int ndigs(char c);
int read_numbers(const char *str, int n, char *dest);
void init_memory(const char *fn, Memory *m);
void show_time(const char *txt);
void output(double *r, unsigned long nr, const char *filename);

#endif
