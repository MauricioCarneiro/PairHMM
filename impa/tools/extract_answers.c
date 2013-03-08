#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_HAPLOTYPE_LEN 10240
#define MAX_READ_LEN 1024

void parse_line(char *l, char *hap, char *read, char *qual, char *ins, char *del, char *cont, int *prev, char *bool_value, double *ans)
{
	memset(hap, 0, MAX_HAPLOTYPE_LEN);
	memset(read, 0, MAX_READ_LEN);
	memset(qual, 0, MAX_READ_LEN);
	memset(ins, 0, MAX_READ_LEN);
	memset(del, 0, MAX_READ_LEN);
	memset(cont, 0, MAX_READ_LEN);
	memset(bool_value, 0, 20);
    sscanf(l, "%s %s %s %s %s %s %d %s %lf\n", hap, read, qual, ins, del, cont, prev, bool_value, ans);
}

void convert(const char *orig)
{
	FILE *from;
	char *l = NULL, haplotype[MAX_HAPLOTYPE_LEN], read[MAX_READ_LEN], qual[MAX_READ_LEN], ins[MAX_READ_LEN], del[MAX_READ_LEN], cont[MAX_READ_LEN], bool_value[20];
	size_t sz = 0;
	int prev;
	double ans;

	from = fopen(orig, "r");

    while (getline(&l, &sz, from) != -1)
	{
		parse_line(l, haplotype, read, qual, ins, del, cont, &prev, bool_value, &ans);
		printf("%e\n", ans); 
	}

	fclose(from);
}

void check_args(int argc, char **argv)
{
	(void) argv;
	if (argc != 2)
	{
		printf("Usage: <binary> <input file>\n");
		exit(0);
	}
}

int main(int argc, char **argv)
{
    check_args(argc, argv);
	convert(argv[1]);
    return 0;
}

