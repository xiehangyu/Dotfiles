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
void interpolationexpression(vector<int> & xs, vector<double> & fs, vector<double> &M, string filename)
{ //calculate the coefficients of the cubic spline, the information is stored in the file
    double a0,a1,a2,a3;
    ofstream fp;//store the coefficients of the cubic spline for later drawing.
    vector<double>h;
    int n=xs.size()-1;
    h.resize(n);
    ofstream fp;
    fp.open(filename);
    if(!fp.is_open())
    {
        cout<<"open file failed"<<endl;
        return;
    }
    cout<<"the cubic spline interpolation expression is:"<<endl;
    fp<<"a3\t"<<"a2\t"<<"a1\t"<<"a0"<<endl;
    for(int i=0;i<=n-1;i++)
    {
        h[i]=xs[i+1]-xs[i];
        a3 = (M[i+1]-M[i])/(6*h[i]);
        a2 = (xs[i+1]*M[i]-xs[i]*M[i+1])/(2*h[i]);
        a1 = (3*(xs[i]*xs[i]*M[i+1]-xs[i+1]*xs[i+1]*M[i])+6*(fs[i+1]-fs[i])-h[i]*h[i]*(M[i+1]-M[i]))/(6*h[i]);
        a0 = (xs[i+1]*xs[i+1]*xs[i+1]*M[i]-xs[i]*xs[i]*xs[i]*M[i+1]+6*(xs[i+1]*fs[i]-xs[i]*fs[i+1])-h[i]*h[i]*(xs[i+1]*M[i]-xs[i]*M[i+1]))/(6*h[i]);
        fp<<a3<<"\t"<<a2<<"\t"<<a1<<"\t"<<a0<<endl;
        cout<<"if x \u2208 ["<<xs[i]<<","<<xs[i+1]<<"],S(x)="<<a3<<"x^3+"<<a2<<"x^2+"<<a1<<"x+"<<a0<<endl;
    }
    fp.close();
}
int main()
{
    vector<int> xs;
    vector<double> fs;
    vector<double> M;
    readdata(xs,fs);
    spline(xs,fs,M);
    interpolationexpression(xs,fs,M,"./data/interpolationexpression_origin.txt");
    cout<<"After the change, the result is:"<<endl;
    xs[9]=0;
    fs[9]=10;
    M.clear();
    spline(xs,fs,M);
    interpolationexpression(xs,fs,M,"./data/interpolationexpression_change.txt");
}