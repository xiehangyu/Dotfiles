#include <stdio.h>
#include <iostream>
#include <random>
using namespace std;
double (*fpp(int N))(double x)
{
	return pow(x,N);
}
int main()
{
	int *p;
	p=new int[7]{1,2,3,4,5,6,7};
	cout<<sizeof(p)/(sizeof(p[0])+0.0)<<endl;
	cout<<endl;
	cout<<p[6]<<endl;
	cout<<fpp(4)(2)<<endl;
}
