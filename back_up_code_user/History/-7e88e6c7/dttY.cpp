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
        if(k==j)
        {
            continue;
        }
        for(int m=1;m<=L-1;m++)
        {
            if(m==k||m==j) 
            {
                continue;
            }
            result += 1/(pow(periodicdistance1D(0,k,L),alpha)*pow(periodicdistance1D(k,m,L),alpha)*pow(periodicdistance1D(j,m,L),alpha));
        }
    }
    return result;
}

#include <iostream>
#include <vector>
#include <cmath>

struct LinearFitResult {
    double m;     // 斜率
    double R2;    // R^2 (决定系数)
};

void LlinearFit(const std::vector<double>& x, const std::vector<double>& y, double & m, double & R2) 
{
    int n = x.size();
    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;

    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }

    m = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    double b = (sum_y - m * sum_x) / n;
    double SST = sum_y2 - sum_y * sum_y / n;
    double SSR = sum_xy - sum_x * sum_y / n;
    double SSE = sum_x2 - sum_x * sum_x / n;
    R2 = 1-SSR / SST;

}

