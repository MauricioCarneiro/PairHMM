#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <string>
#include <sstream>
#include <iostream>

#include "common.h"

double right_now(void)
{
    struct timeval v;
    gettimeofday(&v, (struct timezone *) NULL);
    return v.tv_sec + v.tv_usec/1.0e6;
}

void init_memory(const char *fn, Memory *m)
{
    char *l = NULL;
    int lastptr;
    size_t sz = 0;
    FILE *i;
    unsigned long nG = 0, nR = 0, nH = 0, tr, th, R = 0, H = 0, x = 0;
    unsigned long y = 0, chunk_sz = 0, nres = 0, curr_g = 0, curr_r = 0;
    unsigned long curr_h = 0, curr_res = 0, gr = 0, r = 0, h = 0;

	std::string qual, ins, del, cont;

    i = fopen(fn, "r");
    while (getline(&l, &sz, i) != -1)
    {
        if (strcmp(l, "GROUP\n"))
            break;
        nG++;
        if (getline(&l, &sz, i) == -1)
            return;
        sscanf(l, "%lu %lu\n", &tr, &th);
        nR += tr;
        nH += th;
        nres += (tr * th);
        for (x = 0; x < tr; x++)
        {
            if (getline(&l, &sz, i) == -1) return;
            R = 0; while (l[R] != ' ') R++;
            chunk_sz += (R+1) * 5;
        }
        for (x = 0; x < th; x++)
        {
            if (getline(&l, &sz, i) == -1) return;
            H = 0; while (l[H] != ' ' && l[H] != '\n') H++;
            chunk_sz += (H+1);
        }
    }

    fclose(i);

    m->ng = nG;
    m->g = (Group *) malloc(m->ng * sizeof(Group));

    m->nh = nH;
    m->h = (Haplotype *) malloc(m->nh * sizeof(Haplotype));

    m->nr = nR;
    m->r = (ReadSequence *) malloc(m->nr * sizeof(ReadSequence));

    m->chunk_sz = chunk_sz;
    m->chunk = (char *) malloc(chunk_sz);
    bzero(m->chunk, chunk_sz);

    m->nres = nres;
    m->res = (BIGGEST_NUMBER_REPRESENTATION *) malloc(m->nres * sizeof(BIGGEST_NUMBER_REPRESENTATION));
    m->cmpH = (unsigned long *) malloc(m->nres * sizeof(unsigned long));
    m->cmpR = (unsigned long *) malloc(m->nres * sizeof(unsigned long));

    i = fopen(fn, "r");
    nG = 0;
    curr_g = 0; curr_r = 0; curr_h = 0;
    lastptr = 0;
    while (getline(&l, &sz, i) != -1)
    {
        if (strcmp(l, "GROUP\n")) break;
        nG++;
        if (getline(&l, &sz, i) == -1) return;

        curr_g = nG - 1;
        sscanf(l, "%lu %lu\n", &(m->g[curr_g].nR), &(m->g[curr_g].nH));

        m->g[curr_g].fstR = curr_r;
        m->g[curr_g].fstH = curr_h;

        for (x = 0; x < m->g[curr_g].nR; x++)
        {
            if (getline(&l, &sz, i) == -1) return;

			// read the readsequence curr_r

            R = 0; while (l[R] != ' ') R++;
			// R is the number of characters in the read sequence
            m->r[curr_r].R = R; 
            m->r[curr_r].r = lastptr; // m->r[curr_r].r is the position in the chunk in which the readsequence starts. In the chunk, the readsequence has R characters followed by a '\0' character.
            m->r[curr_r].qidc = lastptr + (R+1); // m->r[curr_r].qidc is the position in the chunk in which the readsequence starts
            sscanf(l, "%s ", m->chunk + lastptr); // read r. Put the R characters into the chunk, starting at m->r[curr_r].r
            y = R+1;

			// read qual, ins, del, cont and put them in an interleaved way into the chunk.
			std::istringstream iss(l + y);
			iss >> qual >> ins >> del >> cont;
			std::istringstream issq(qual);
			std::istringstream issi(ins);
			std::istringstream issd(del);
			std::istringstream issc(cont);
			for (int j = 0; j < R; j++)
			{
				if (issq.peek() == ',') issq.ignore();
				if (issi.peek() == ',') issi.ignore();
				if (issd.peek() == ',') issd.ignore();
				if (issc.peek() == ',') issc.ignore();
				int intq, inti, intd, intc;
				issq >> intq; if (intq < 6) intq = 6;
				issi >> inti;
				issd >> intd;
				issc >> intc;
				m->chunk[lastptr + (R+1) + 4 * j + 0] = (char) intq;
				m->chunk[lastptr + (R+1) + 4 * j + 1] = (char) inti;
				m->chunk[lastptr + (R+1) + 4 * j + 2] = (char) intd;
				m->chunk[lastptr + (R+1) + 4 * j + 3] = (char) intc;
			}

            lastptr += (R+1) * 5;
            curr_r++;
        }
        for (x = 0; x < m->g[curr_g].nH; x++)
        {
            if (getline(&l, &sz, i) == -1) return;
            H = 0; while (l[H] != ' ' && l[H] != '\n') H++;

            m->h[curr_h].H = H;
            m->h[curr_h].h = lastptr;
            sscanf(l, "%s\n", m->chunk + lastptr);

            lastptr += H+1;
            curr_h++;
        }
    }

    fclose(i);

    curr_res = 0;
    for (gr = 0; gr < m->ng; gr++)
        for (r = 0; r < m->g[gr].nR; r++)
            for (h = 0; h < m->g[gr].nH; h++)
            {
                m->cmpH[curr_res] = m->g[gr].fstH + h;
                m->cmpR[curr_res] = m->g[gr].fstR + r;
                curr_res++;
            }
    return;
}

void output(BIGGEST_NUMBER_REPRESENTATION *r, unsigned long nr, const char *filename)
{
    FILE *f;
    unsigned long x;
    f = fopen(filename, "w");
    for (x = 0; x < nr; x++)
        fprintf(f, "%.18E\n", r[x]);
    fclose(f);
}

