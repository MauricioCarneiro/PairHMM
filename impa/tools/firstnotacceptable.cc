#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using std::ifstream;
using std::ofstream;
using std::string;
using std::istringstream;

typedef const char cch;
typedef double dbl;

void findfst(cch *fn1, cch *fn2, int &vpos, dbl &v, dbl &v1, dbl &v2, dbl bound)
{
	ifstream f1(fn1, ifstream::in);
	ifstream f2(fn2, ifstream::in);
	string s1, s2;
	bool found = false;
	vpos = 0;

	while (getline(f1, s1) && getline(f2, s2))
	{
		vpos++;
        istringstream l1(s1);
        istringstream l2(s2);
		l1 >> v1;
		l2 >> v2;
		v = fabs(v1 - v2);
		if (v > bound)
		{
			found = true;
			break;
		}
	}

	if (!found)
	{
		v1 = 0;
		v2 = 0;
		v = 0.0;
		vpos = -1;
	}

	f1.close();
	f2.close();
}

void check_args(int argc, char **argv)
{
    FILE *f1, *f2, *f3;

    if (argc != 6)
    {
        printf("Usage: firstnotacceptable <file1> <file2> <highest acceptable difference> <txt input> <txt output>\n");
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

    if (!(f3 = fopen(argv[4], "r")))
    {
        printf("Cannot open %s\n", argv[4]);
        exit(0);
    }
    fclose(f3);

    return;
}

void print_isolated_case(const char *fin, int row, const char *fout)
{
	ifstream f1(fin, ifstream::in);
	string s1;

	for (int x = 1; x < row; x++)
		getline(f1, s1);
	getline(f1, s1);
	f1.close();

	ofstream f2(fout);
	f2 << s1;
	f2.close();
}


int main(int argc, char **argv)
{
    int row;
    double diff, v1, v2, bound;
    check_args(argc, argv);
	bound = atof(argv[3]);
    printf("Looking for the first p2p difference higher than %e\n", bound);
    printf("--------------------------------------------------------\n");
    findfst(argv[1], argv[2], row, diff, v1, v2, bound);
    printf("The difference is: %e\n", diff);
    printf("Occurs at row    : %d\n", row);
    printf("%s value: %e\n", argv[1], v1);
    printf("%s value: %e\n", argv[2], v2);
    printf("--------------------------------------------------------\n");
    printf("Isolated case printed to %s\n", argv[5]);
	print_isolated_case(argv[4], row, argv[5]);
}

