#include <string>
#include <iostream>
#include <vector>

using std::string;
using std::vector;
using std::cin;
using std::cout;

int main()
{
	int num_of_read_sequences, num_of_haplotypes;
	vector<string> haps, rss, qs, is, ds, cs;
	string hap, rs, q, i, d, c;

	while ((cin >> num_of_read_sequences >> num_of_haplotypes).good())
	{
		haps.clear();
		rss.clear();
		qs.clear();
		is.clear();
		ds.clear();
		cs.clear();

		for (int r = 0; r < num_of_read_sequences; r++)
		{
			cin >> rs >> q >> i >> d >> c;
			rss.push_back(rs);
			qs.push_back(q);
			is.push_back(i);
			ds.push_back(d);
			cs.push_back(c);
		}

		for (int h = 0; h < num_of_haplotypes; h++)
		{
			cin >> hap;
			haps.push_back(hap);
		}

		for (int h = 0; h < num_of_haplotypes; h++)
			for (int r = 0; r < num_of_read_sequences; r++)
				cout << haps[h] << " " << rss[r] << " " << qs[r] << " " << is[r] << " " << ds[r] << " " << cs[r] << " 0 true\n";

	}
}
