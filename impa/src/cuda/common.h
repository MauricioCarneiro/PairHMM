#ifndef COMMON_H
#define COMMON_H

#define BIGGEST_NUMBER_REPRESENTATION double

typedef struct
{
	unsigned long nR, nH, fstR, fstH;
} Group;

typedef struct
{
	unsigned long R;
	int r, qidc;
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
    BIGGEST_NUMBER_REPRESENTATION *res;
	char *flag;
    unsigned long ng, nh, nr, chunk_sz, nres;
} Memory;

double right_now(void);
void init_memory(const char *fn, Memory *m);
void output(BIGGEST_NUMBER_REPRESENTATION *r, unsigned long nr, const char *filename);

#endif
