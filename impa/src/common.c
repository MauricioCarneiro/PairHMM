#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "common.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#include <sys/time.h>
#endif

double right_now(void)
{
#ifdef _WIN32
    LARGE_INTEGER counter, freq;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&freq);
    return (1.0*counter.QuadPart)/(1.0*freq.QuadPart);
#else
    struct timeval v;
    gettimeofday(&v, (struct timezone *) NULL);
    return v.tv_sec + v.tv_usec/1.0e6;
#endif
}

int ndigs(char c)
{
    int t = 0;
    if (c < 0) { c *= -1; t++; }
    if ((0 <= c) && (c < 10)) return (t + 1);
    if ((10 <= c) && (c < 100)) return (t + 2);
    return (t + 3);
}

int read_numbers(const char *str, int n, char *dest)
{
    int pos = 0, x, temp;
    sscanf(str, "%d", &temp);
    dest[0] = temp;
    pos += ndigs(*dest);
    for (x = 1; x < n; x++)
    {
        sscanf(str + pos, ",%d", &temp);
        dest[x] = (char)temp;
		dest[x] = ((dest[x]) & 127);
        pos += ndigs(dest[x]) + 1;
    }
    return pos;
}

void init_memory(const char *fn, Memory *m)
{
    char *l = NULL;
    int lastptr;
    size_t sz = 0;
    FILE *i;
    unsigned long nG = 0, nR = 0, nH = 0, tr, th, R = 0, H = 0, x = 0;
    unsigned long y = 0, t = 0, chunk_sz = 0, nres = 0, curr_g = 0, curr_r = 0;
    unsigned long curr_h = 0, curr_res = 0, gr = 0, r = 0, h = 0;

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
	m->res = (double *) malloc(m->nres * sizeof(double));
	m->cmpH = (unsigned long *) malloc(m->nres * sizeof(unsigned long));
	m->cmpR = (unsigned long *) malloc(m->nres * sizeof(unsigned long));

    m->phred_to_prob = (double *)malloc(256 * sizeof(double));
    for (t = 0; t < 256; t++)
        m->phred_to_prob[t] = pow(10.0, -((double) t) / 10.0);

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
            R = 0; while (l[R] != ' ') R++;

            m->r[curr_r].R = R;
            m->r[curr_r].r = lastptr;
            m->r[curr_r].qual = lastptr + (R+1);
            m->r[curr_r].ins = lastptr + (R+1) * 2;
            m->r[curr_r].del = lastptr + (R+1) * 3;
            m->r[curr_r].cont = lastptr + (R+1) * 4;
            sscanf(l, "%s ", m->chunk + lastptr);
            y = R+1;
            y += read_numbers(l + y, R, m->chunk + lastptr + (R+1)); y++;
            for (t = 0; t < R; t++) 
                m->chunk[m->r[curr_r].qual + t] = (m->chunk[m->r[curr_r].qual + t] < 6) ? 6 : m->chunk[m->r[curr_r].qual + t];
            y += read_numbers(l + y, R, m->chunk + lastptr + (R+1)*2); y++;
            y += read_numbers(l + y, R, m->chunk + lastptr + (R+1)*3); y++;
            read_numbers(l + y, R, m->chunk + lastptr + (R+1)*4);

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

void show_time(const char *txt)
{
    static int N = 0;
    static double previous;
    static double start;
    double now;

    if (N == 0)
    {
        start = right_now();
        previous = start;
    }

    now = right_now();

    printf("Prev: %e sec   Tot: %e sec ===   %s\n", now - previous, now - start, txt);

    N++;
    previous = now;
    return;
}

void output(double *r, unsigned long nr, const char *filename)
{
    FILE *f;
    unsigned long x;
    f = fopen(filename, "w");
    for (x = 0; x < nr; x++)
        fprintf(f, "%.18E\n", r[x]);
    fclose(f);
}


