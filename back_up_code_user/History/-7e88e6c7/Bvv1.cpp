#include <iostream>
#include <cmath>
/*
guess when alpha \leq 2/3, \sim L^(2-3alpha)
when 2/3<alpha \leq 1, \sim |ri-rj|^(2-3alpha)
when alpha > 1 \sim |ri-rj|^(-alpha) 
*/
using namespace std;

double periodicdistance1D(int x1, int x2, int L)
{
    auto d = abs(x1-x2);
    return min(d, L-d);
}

double D1casecalculation(int j, double alpha, int L)
{
    double result=0;
    for(int k=1;k<=L-1;k++)
    {
        for(int m=1;m<=L-1;m++)
        {
            
        }
    }
}