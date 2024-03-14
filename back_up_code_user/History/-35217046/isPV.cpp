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
int main()
{
	double a[100];
	double b[100];
	double c[100];
	for(int i=0;i<=99;i++)
	{
		a[i]=i;
		b[i]=i;
		c[i]=a[i]+b[i];
	}
	printf("success!");
}
