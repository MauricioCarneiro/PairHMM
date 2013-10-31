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
    BIGGEST_NUMBER_REPRESENTATION *res;
	char *flag;
    BIGGEST_NUMBER_REPRESENTATION *phred_to_prob; /* not used */
    unsigned long ng, nh, nr, chunk_sz, nres;
} Memory;

double right_now(void);
int ndigs(char c);
int read_numbers(const char *str, int n, char *dest);
void init_memory(const char *fn, Memory *m);
void show_time(const char *txt);
void output(BIGGEST_NUMBER_REPRESENTATION *r, unsigned long nr, const char *filename);

#endif
