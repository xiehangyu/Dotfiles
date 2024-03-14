#include <stdio.h>
#include <iostream>
#include <random>
#include <eigen3/Eigen/Core>
using namespace std;
class mypower
{
	public:
		mypower(int n):n(n){};
		double operator()(double x)
{
	return pow(x,n);
}
	private:
		int n;
};
void stencil_softpipe(double *Y, double alpha, const int N)
{
    double new3,new2,new1,new4,new5,new6;
    new1=Y[1];
    new2=Y[2];
    int upperbound=(N-2)/4*4;
    for(int i=1;i<N-2;i+=4)
    {
        new3=Y[i+2];
        Y[i]=alpha * (new2+Y[i-1])+new1;
        new2=new3;
        new1=new2;
    }
    Y[N-2]=alpha * (new2+Y[N-3])+new1;
}
int main()
{
	double Y[5]={1,2,3,4,5};
	stencil_softpipe(Y,1,5);
	for(int i=0;i<5;i++)
	{
		cout<<Y[i]<<endl;
	}
}
