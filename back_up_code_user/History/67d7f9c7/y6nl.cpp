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
    vector<double> b;
    int gapnumber=n-1;
    h.resize(gapnumber);
    b.resize(gapnumber);
    for(int i=0;i<=gapnumber-1;i++)
    {
        h[i]=xs[i+1]-xs[i];
        b[i]=6.0/(h[i])*(fs[i+1]-fs[i]);
    }
    //using LU decomposition to solve the linear equation
    vector<double> u;
    vector<double> v;
    u.resize(gapnumber-1);
    v.resize(gapnumber-1);
    u[0]=2.0*(h[0]+h[1]);
    v[0]=b[1]-b[0];
    for(int i=1;i<=gapnumber-2;i++)
    {
        u[i]=2.0*(h[i]+h[i+1])-h[i]*h[i]/u[i-1];
        v[i]=b[i+1]-b[i]-h[i]*v[i-1]/u[i-1];
    }
    for(int i=gapnumber-1;i>=1;i--)
    {
        M[i]=(v[i-1]-h[i]*M[i+1])/u[i-1];
    }
}
void interpolationexpression(vector<int> & xs, vector<double> & fs, vector<double> &M)
{
    ofstream fp;//store the coefficients of the cubic spline for later drawing.
    vector<double>h;
    int n=xs.size()-1;
    h.resize(n);
    fp.open("./data/interpolationexpression.txt");
    if(!fp.is_open())
    {
        cout<<"open file failed"<<endl;
        return;
    }
    for(int i=0;i<=n-1;i++)
    {
        h[i]=xs[i+1]-xs[i];
    }
}