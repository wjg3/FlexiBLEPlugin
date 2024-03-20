#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

int main()
{
    fstream fin("water_coor.txt", ios::in);
    fstream fout("qm.txt", ios::out);
    vector<double> r;
    for (int i = 0; i < 600; i++)
    {
        double x, y, z;
        fin >> x >> y >> z;
        if (i % 3 == 0)
        {
            double rNow = sqrt(x * x + y * y + z * z);
            r.emplace_back(rNow);
        }
    }
    vector<double> r_sorted = r;
    stable_sort(r_sorted.begin(), r_sorted.end());
    for (int i = 0; i < 200; i++)
    {
        if (r[i] <= r_sorted[19])
        {
            fout << i * 3 << ",";
            fout << i * 3 + 1 << ",";
            fout << i * 3 + 2 << ",";
        }
    }
}