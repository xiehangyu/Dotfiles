#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
using namespace std;
void readdata(vector<int> & xs, vector<double> &fs)
{
    //read data from the file
    ifstream fp("./data/point.txt");
    if (!fp.is_open())
    {
        cout << "open file failed" << endl;
        return;
    }
    int x;
    double f;
    while (fp >> x >> f)
    {
        xs.push_back(x);
        fs.push_back(f);
    }
}
void spline(vector<int> & xs, vector<double> & fs, vector<double> & M)
{
    //the cubic spline function, using natural condition
    int n=xs.size();
    M.resize(n);
    M[0]=0;
    M[n-1]=0;//natural condition
    vector<double> h;
    h.resize(n-1);
    for(int i=0;i<=n-2;i++)
}