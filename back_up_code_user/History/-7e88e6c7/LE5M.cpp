#include <iostream>
#include <cmath>

using namespace std;

double periodicdistance1D(int x1, int x2, int L)
{
    int dx = x1 - x2;
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