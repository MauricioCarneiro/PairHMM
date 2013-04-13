#include <iostream>
#include <cstdlib>
#include <ctime>

using std::cout;

#define MB 900

int main(void)
{
	char *d_c, *h_c;
	srand(time(NULL));
	cudaMalloc(&d_c, MB * 1024 * 1024);
	h_c = (char *)malloc(MB * 1024 * 1024);

	cout << "d_c = ***" << (void *)d_c << "***\n";
	cout << "h_c = ***" << (void *)h_c << "***\n";

	if (d_c == 0)
	{
		cout << "cudaMalloc error\n";
		exit(0); 
	}

	for (int x = 0; x < MB * 1024 * 1024; x++)
	{
		h_c[x] = rand() % 90 + 32;
		if (x % 1000 == 0)
			cout << x << "\n";
	}
	cudaMemcpy(d_c, h_c, MB * 1024 * 1024, cudaMemcpyHostToDevice);
	cudaFree(d_c);
	free(h_c);
}
