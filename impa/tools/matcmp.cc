#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

void compute_biggest_difference(const char *fn1, const char *fn2, int h, int w, int *ret_row, int *ret_col, float *ret, float *n1, float *n2)
{
	std::ifstream if1(fn1 , std::ifstream::in);
	std::ifstream if2(fn2 , std::ifstream::in);
	std::string s1, s2;

    float f1, f2, diff;
	int r, c;

    *ret = -1;

    for (r = 0; r < h; r++)
        for (c = 0; c < w; c++)
		{
			if1 >> s1;
			if2 >> s2;
			
			f1 = strtof(s1.c_str(), NULL);
			f2 = strtof(s2.c_str(), NULL);

            diff = fabs(f1 - f2);
			if (fabs(f1) == INFINITY && fabs(f2) == INFINITY)
				diff = f1 * f2 > 0 ? 0 : -INFINITY;
			
            if (diff > *ret)
            {
                *ret = diff;
                *ret_row = r;
                *ret_col = c;
                *n1 = f1;
                *n2 = f2;
            }
		}

    return;
}

int get_dimensions(const char *fn, int *height, int *width)
{
    FILE *f;
    int i, tc;
    char buf[10240];
	std::string s1, s2;

    *height = 0;
    *width = 0;

    f = fopen(fn, "r");

    for (i = fscanf(f, "%[^\n]", buf); i != EOF; i = fscanf(f, "\n%[^\n]", buf))
    {
        std::istringstream line(buf);
        tc = 0;
        while (line >> s1)
            tc++;
        
        if (*height == 0)
            (*width) = tc;
        else
            if (tc != *width)
            {
                printf("%s: line %d has %d columns while previous ones had %d columns.\n", fn, (*height)+1, tc, *width);
                exit(0);
            }

        (*height)++;
    }

    fclose(f);
    return 0;
}

void check_args(int argc, char **argv)
{
    FILE *f1, *f2;

    if (argc != 3)
    {
        printf("Usage: compare <file1> <file2>\n");
        exit(0);
    }

    if (!(f1 = fopen(argv[1], "r")))
    {
        printf("Cannot open %s\n", argv[1]);
        exit(0);
    }
    fclose(f1);

    if (!(f2 = fopen(argv[2], "r")))
    {
        printf("Cannot open %s\n", argv[2]);
        exit(0);
    }
    fclose(f2);
    return;
}

int main(int argc, char **argv)
{
    int h1, w1, h2, w2, maxdif_r, maxdif_c;
    float maxdif, n1, n2;

    check_args(argc, argv);
    printf("Getting dimensions of %s: ", argv[1]);
    get_dimensions(argv[1], &h1, &w1);
    printf("Done! H = %d; W = %d\n", h1, w1);
    printf("Getting dimensions of %s: ", argv[2]);
    get_dimensions(argv[2], &h2, &w2);
    printf("Done! H = %d; W = %d\n", h2, w2);
    if (h1 != h2 || w1 != w2)
    {
        printf("Dimensions of both files are different!!!\n");
        exit(0);
    }
    printf("Looking for the biggest point-to-point difference\n");
    printf("--------------------------------------------------------\n");
    compute_biggest_difference(argv[1], argv[2], h1, w1, &maxdif_r, &maxdif_c, &maxdif, &n1, &n2);
    printf("The biggest difference is: %f\n", maxdif);
    printf("Occurs at:                 %d, %d\n", maxdif_r, maxdif_c);
    printf("%s value is:        %f\n", argv[1], n1);
    printf("%s value is:        %f\n", argv[2], n2);
    printf("--------------------------------------------------------\n");
}

