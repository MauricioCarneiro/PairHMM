#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "common.h"

void statistics(Memory *m)
{
	unsigned long g, h, r, R, H, mem, maxmem = 0, maxmemR, maxmemH, maxR = 0, maxH = 0;
	unsigned long maxcomps = 0, maxcompsnh, maxcompsnr;
	for (g = 0; g < m->ng; g++)
	{
		if (m->g[g].nH * m->g[g].nR > maxcomps)
		{
			maxcomps = m->g[g].nH * m->g[g].nR;
			maxcompsnh = m->g[g].nH;
			maxcompsnr = m->g[g].nR;
		}
		for (h = m->g[g].fstH; h < m->g[g].fstH + m->g[g].nH - 1; h++)
			for (r = m->g[g].fstR; r < m->g[g].fstR + m->g[g].nR - 1; r++)
			{
				R = m->r[r].R;
				H = m->h[h].H;
				mem = ((R+1) * 15 + 256) * 8 + (H+1) + (R+1) * 5;
				if (mem > maxmem)
				{
					maxmem = mem;
					maxmemR = R;
					maxmemH = H;
				}
				if (R > maxR) maxR = R;
				if (H > maxH) maxH = H;
			}
	}
	printf("The maximum amount of memory needed by a comparison is: %lu\n", maxmem);
	printf("   R: %lu\n", maxmemR);
	printf("   H: %lu\n", maxmemH);
	printf("Maximum R: %lu\n", maxR);
	printf("Maximum H: %lu\n", maxH);
	printf("Maximum number of comparisons in a group: %lu\n", maxcomps);
	printf("   nR: %lu\n", maxcompsnr);
	printf("   nH: %lu\n", maxcompsnh);
}

int main(int argc, char **argv)
{
    Memory m;

    if (argc != 2)
    {
        printf("\nUsage: <binary> <input file>\n\n");
        exit(0);
    }

    init_memory(argv[1], &m);
	statistics(&m);
    return 0;
}

