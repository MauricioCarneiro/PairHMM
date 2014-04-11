#ifndef TIMING_H
#define TIMING_H

#include <ctime>
#include <sys/time.h>
#include <string>

using std::string;

class Timing
{
public:
	Timing(string t) : st(now()), tot(0.0), title(t)
	{
	}

	~Timing()
	{
		std::cerr << title << tot << " seconds\n";
	}

	void start()
	{
		st = now();
	}

	void acc()
	{
		tot += (now() - st);
	}

private:
	static double now()
	{
		struct timeval v;
		gettimeofday(&v, (struct timezone *) NULL);
		return v.tv_sec + v.tv_usec/1.0e6;
	}

	double st, tot;
	string title;
};

#endif
