#include <iostream>
#include <cmath>

using namespace std;

double periodicdistance1D(double x1, double x2, int L)
{
    double dx = x1 - x2;
    if (dx > L/2)
    {
        dx -= L;
    }
    else if (dx < -L/2)
    {
        dx += L;
    }
    return dx;
}

double 1Dcase(int i, int j, double alpha, int L)