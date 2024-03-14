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
	double sum=0;
	for(int i=0;i<=1000000000;i++)
	{
		sum+=i;
	}
	printf("%f\n",sum);
}
