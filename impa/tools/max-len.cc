#include <string>
#include <iostream>

using std::string;
using std::cin;
using std::cout;

void max_len(int &haplen, int &rslen)
{
	string hap, rs, q, i, d, c, i1, i2;
	haplen = -1;
	rslen = -1;

	while ((cin >> hap >> rs >> q >> i >> d >> c >> i1 >> i2).good())
	{
		int szhap = hap.size();
		int szrs = rs.size();

		if (szhap > haplen)
			haplen = szhap;

		if (szrs > rslen)
			rslen = szrs;

		string("").swap(hap);
		string("").swap(rs);
		string("").swap(q);
		string("").swap(i);
		string("").swap(d);
		string("").swap(c);
		string("").swap(i1);
		string("").swap(i2);
	}
}

int main(int argc, char *argv[])
{
	int maxhap, maxrs;
	max_len(maxhap, maxrs);
	cout << "max rslen: " << maxrs << "\n";
	cout << "max haplen: " << maxhap << "\n";
		
	return 0;
}

