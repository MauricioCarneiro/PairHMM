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
using std::ostringstream;
using std::cout;

#define DELIM ','

string hms2secs(string cad)
{
	if (cad == "")
		return string("");

	const char *hms = cad.c_str();
	ostringstream curr("");
	double ret = 0;

	for (int x = 0; x < (int)cad.size(); x++)
	{
		if (hms[x] == ' ')
			continue;

		if (hms[x] == 'h')
		{
			ret += (atoi(curr.str().c_str()) * 3600);
			curr.str("");
			curr.clear();
			continue;
		}
		if (hms[x] == 'm')
		{
			ret += (atoi(curr.str().c_str()) * 60);
			curr.str("");
			curr.clear();
			continue;
		}
		if (hms[x] == 's')
		{
			ret += (atof(curr.str().c_str()));
			curr.str("");
			curr.clear();
			continue;
		}
		if ((hms[x] >= '0' && hms[x] <= '9') || (hms[x] == '.'))
		{
			curr << hms[x];
		}
	}
	curr.str("");
	curr.clear();
	curr << ret;
	return curr.str();
}

int main(int argc, char **argv)
{
	(void) argc;
	(void) argv;

	string line;
	//YEAR;MONTH;DAY;HOUR;MINUTE;LANGUAGE;TEST;HOST;GPU_MODEL;SIZE;COMMIT_ID;INIT_MEMORY;MALLOCS;COMPUTATION;KERNEL;OUTPUT;TOT_(measured_by_program);TOT_(measured_by_time);BIGGEST_ERROR
	string year, month, day, hour, minute, language, test, host, gpu_model, size, commit_id;
	string init_memory, mallocs, computation, kernel, output, tot_program, tot_time, error;

	int intyear, intmonth, intday, inthour, intminute, intsecond=0;
	char date[15];

	int linenum=0;
	while (getline(std::cin, line))
	{
		istringstream iss(line);
		getline(iss, year, ';');
		getline(iss, month, ';');
		getline(iss, day, ';');
		getline(iss, hour, ';');
		getline(iss, minute, ';');
		getline(iss, language, ';');
		getline(iss, test, ';');
		getline(iss, host, ';');
		getline(iss, gpu_model, ';');
		getline(iss, size, ';');
		getline(iss, commit_id, ';');
		getline(iss, init_memory, ';');
		getline(iss, mallocs, ';');
		getline(iss, computation, ';');
		getline(iss, kernel, ';');
		getline(iss, output, ';');
		getline(iss, tot_program, ';');
		getline(iss, tot_time, ';');
		getline(iss, error, ';');

		if (linenum++ == 0)
		{
			cout << "DATE" << DELIM;
			cout << "LANGUAGE" << DELIM;
			cout << "TEST" << DELIM;
			cout << "HOST" << DELIM;
			cout << "GPU_MODEL" << DELIM;
			cout << "SIZE" << DELIM;
			cout << "COMMIT_ID" << DELIM;
			cout << "INIT_MEMORY" << DELIM;
			cout << "MALLOCS" << DELIM;
			cout << "COMPUTATION" << DELIM;
			cout << "KERNEL" << DELIM;
			cout << "OUTPUT" << DELIM;
			cout << "TOT_(measured_by_program)" << DELIM;
			cout << "TOT_(measured_by_time)" << DELIM;
			cout << "BIGGEST_ERROR\n";
			continue;
		}

		date[14]=0;
		intyear = atoi(year.c_str());
		intmonth = atoi(month.c_str());
		intday = atoi(day.c_str());
		inthour = atoi(hour.c_str());
		intminute = atoi(minute.c_str());
		sprintf(date, "%04d%02d%02d%02d%02d%02d", intyear, intmonth, intday, inthour, intminute, intsecond);

		cout << date << DELIM;
		cout << language << DELIM;
		cout << test << DELIM;
		cout << host << DELIM;
		cout << gpu_model << DELIM;
		cout << size << DELIM;
		cout << commit_id << DELIM;
		cout << hms2secs(init_memory) << DELIM;
		cout << hms2secs(mallocs) << DELIM;
		cout << hms2secs(computation) << DELIM;
		cout << hms2secs(kernel) << DELIM;
		cout << hms2secs(output) << DELIM;
		cout << hms2secs(tot_program) << DELIM;
		cout << hms2secs(tot_time) << DELIM;
		cout << error << "\n";
	}
}

