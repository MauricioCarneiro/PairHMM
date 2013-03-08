#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	int h, m;
	double s, t;
	(void) argc;

	if (argc != 3)
		return 0;

	if ((strcmp(argv[1], "ms")) && (strcmp(argv[1], "s")))
		return 0;

	t = atof(argv[2]);

	if (!strcmp(argv[1], "ms"))
		t /= 1000.0;

	h = t / 3600;
	m = (t - h * 3600) / 60;
	s = t - h * 3600 - m * 60;

	if (h > 0) printf("%dh ", h);
	if ((h > 0) || (m>0)) printf("%dm ", m);
	printf("%.3fs", s);
	return 0;
}
