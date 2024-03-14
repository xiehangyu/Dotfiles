#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
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
    #pragma omp parallel for reduction(+:result)
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

void LlinearFit(const std::vector<double>& x, const std::vector<double>& y, double & m, double & b, double & R2) 
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
    b = (sum_y - m * sum_x) / n;
    double SST = sum_y2 - sum_y * sum_y / n;
    std::vector<double> y_hat(n);
    for(int i = 0; i < n; i++)
    {
        y_hat[i] = m * x[i] + b;
    }
    double SSR = 0;
    for(int i = 0; i < n; i++)
    {
        SSR += (y_hat[i] - y[i]) * (y_hat[i] - y[i]);
    }
    R2=1-SSR/SST;

}

void fitting_ij()
{
    int L;
    cout<<"Please input alpha: ";
    double alpha;
    cin>>alpha;
    cout<<"Please input L: ";
    cin>>L;
    vector<double> x;
    vector<double> y;
    for(int j=1;j<=800;j++)
    {
        x.push_back(log(j));
        y.push_back(log(D1casecalculation(j,alpha,L)));
    }
    double m, R2;
    LlinearFit(x,y,m,R2);
    cout<<"m: "<<m<<endl;
    cout<<"R2: "<<R2<<endl;
}

void fitting_L()
{
    cout<<"Please input alpha: ";
    double alpha;
    cin>>alpha;
    cout<<"Please input j: ";
    int j;
    cin>>j;
    vector<double> x;
    vector<double> y;
    for(int L=2000;L<=8000;L++)
    {
        x.push_back(log(L));
        y.push_back(log(D1casecalculation(j,alpha,L)));
    }
    double m, R2;
    LlinearFit(x,y,m,R2);
    cout<<"m: "<<m<<endl;
    cout<<"R2: "<<R2<<endl;
}

int main()
{
    omp_set_num_threads(8);
    char choice;
    do
    {
        cout<<"1 for fitting_ij, 2 for fitting_L, 3 for exit: ";
        cin>>choice;
        switch(choice)
        {
            case '1':
                fitting_ij();
                break;
            case '2':
                fitting_L();
                break;
            case '3':
                break;
            default:
                cout<<"Wrong input!"<<endl;
                break;
        }
    }while(choice!='3');
}