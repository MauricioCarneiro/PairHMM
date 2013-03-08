#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

typedef pair<string, vector<vector<int> > > READ;
typedef string HAPLOTYPE;

struct Iteration
{
    READ r;
    vector<HAPLOTYPE> hts;
};

struct Group
{
    vector<READ> rs;
    vector<HAPLOTYPE> hts;
};

void shownumbers(vector<int> &v)
{
    cout << v[0];
    for (int x = 1; x < (int)v.size(); x++)
        cout << "," << v[x];
}

int normalize(char c)
{
	return ((int) (c - 33));
}

void parseline(string l, HAPLOTYPE &ht, READ &r)
{
    vector<int> vi;
    string qual, ins, del, cont;

    r.second.clear();
    istringstream iss(l);
    iss >> ht >> r.first >> qual >> ins >> del >> cont;

    vi.clear(); 
    for (int x = 0; x < (int)r.first.size(); x++) 
		vi.push_back(normalize((qual.c_str())[x]));
    r.second.push_back(vi);

    vi.clear(); 
    for (int x = 0; x < (int)r.first.size(); x++) 
		vi.push_back(normalize((ins.c_str())[x]));
    r.second.push_back(vi);

    vi.clear(); 
    for (int x = 0; x < (int)r.first.size(); x++) 
		vi.push_back(normalize((del.c_str())[x]));
    r.second.push_back(vi);

    vi.clear(); 
    for (int x = 0; x < (int)r.first.size(); x++) 
		vi.push_back(normalize((cont.c_str())[x]));
    r.second.push_back(vi);

}

void showgroup(Group &group)
{
    cout << "GROUP\n";
    cout << group.rs.size() << " " << group.hts.size() << "\n";
    for (int w = 0; w < (int)group.rs.size(); w++)
    {
        cout << group.rs[w].first << " ";
        shownumbers(group.rs[w].second[0]); cout << " ";
        shownumbers(group.rs[w].second[1]); cout << " ";
        shownumbers(group.rs[w].second[2]); cout << " ";
        shownumbers(group.rs[w].second[3]); cout << "\n";
    }

    for (int w = 0; w < (int)group.hts.size(); w++)
        cout << group.hts[w] << "\n";
}

void flush_group(Group &g)
{
    if (g.rs.size() > 0)
        showgroup(g);
    g.hts.clear();
    g.rs.clear();
}

void process_iteration(Group &g, Iteration &it)
{
    if ((g.hts.size() == 0) || it.hts == g.hts) 
    {
        if (g.hts.size() == 0) 
            g.hts = it.hts;
        g.rs.push_back(it.r);
        return;
    }
    flush_group(g);
    g.hts = it.hts;
    g.rs.push_back(it.r);
}

int main(int argc, char **argv)
{
    ifstream i;
    string l;
    HAPLOTYPE ht;
    READ r;
    Iteration it;
    Group g;

    if (argc != 2) { cout << "\n\nUsage <binary> <input>\n\n"; }

    i.open(argv[1]);
    while (getline(i, l))
    {
        parseline(l, ht, r);

        if ((it.hts.size() == 0) || r == it.r) // EMPTY(it) || SAME_ITER(r, it)
        {
            if (it.hts.size() == 0)
                it.r = r;
            it.hts.push_back(ht);
        }
        else // NOT_EMPTY(it) && NOT_SAME_ITER(r, it)
        {
            process_iteration(g, it);
            it.r = r;
            it.hts.clear();
            it.hts.push_back(ht);
        }
    }
    process_iteration(g, it);
    flush_group(g);
    
    i.close();
    return 0;
}

