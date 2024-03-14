#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>

using namespace itensor;
using namespace std;


static int N;
static SiteSet sites;

double Period_length(int i,int j)
{
    return min(fabs(i-j+0.0,)N-fabs(i-j+0.0));
}

auto LongrangeAntiHeisenberg(double alpha, double h=0.0)
{
    auto ampo=AutoMPO(sites);
    double normfactor=max(pow(N,1-alpha),1.0);
    normfactor=1.0/normfactor;
    for(int i=1;i<=N;i++)
    {
        for(int j=i+1;j<=N;j++)
        {
            ampo+=normfactor/pow(1.0+Period_length(i,j),alpha),"Sx",i,"Sx",j;
            ampo+=normfactor/pow(1.0+Period_length(i,j),alpha),"Sy",i,"Sy",j;
            ampo+=normfactor/pow(1.0+Period_length(i,j),alpha),"Sz",i,"Sz",j;
        }
        ampo += h,"Sz",i;
    }
    auto H=toMPO(ampo,{"Exact=",true});
    printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    return H;
}

auto productstate()
{
    auto state=initState(sites,"Up");
    return MPS(state);
}

