#include <stdio.h>
#include <iostream>
#include <random>
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

	int *p;
	p=new int[7]{1,2,3,4,5,6,7};
	cout<<sizeof(p)/(sizeof(p[0])+0.0)<<endl;
	cout<<endl;
	cout<<p[6]<<endl;
	mypower mp(4);
	mypower mp2=&mp;
	cout<<mp2(2)<<endl;
}
