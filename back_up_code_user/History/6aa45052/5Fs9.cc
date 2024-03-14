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
    return min(fabs(i-j+0.0),N-fabs(i-j+0.0));
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
    auto state=InitState(sites,"Up");
    return MPS(state);
}
auto randomstate()
{
    auto state=InitState(sites,"Up");
    return randomMPS(state);
}
auto correlators2point(int i,int j, string s1, string s2, MPS psi)
{
    auto op_i = op(sites,s1,i);
    auto op_j = op(sites,s2,j);
    psi.position(i);
    auto psidag = dag(psi);
    psidag.prime("Link");
    auto li_1=leftLinkIndex(psi,i);
    auto C=prime(psi(i),li_1)*op_i;
    C *= prime(psidag(i),"Site");
    for(int k=i+1;k<j;k++)
    {
        C *= psi(k);
        C *= psidag(k);
    }
    auto lj=rightLinkIndex(psi,j);
    C*=prime(psi(j),lj)*op_j;
    C*=prime(psidag(j),"Site");
    auto result = eltC(C);
    return result;
}

void imaginarytimeevolutionTDVP(double alpha, double t=0.1, double tend=20, int n_GSE=3, int Korder=3, double h=0)
{
    auto psi=productstate();
    auto H=LongrangeAntiHeisenberg(alpha,h);
    int nsw=round(tend/t);
    auto sweeps=Sweeps(1);
    sweeps.maxdim()=10000;
    sweeps.cutoff()=1E-8;
    sweeps.niter()=50;
    cout<<endl;
    cout<<"The initial energy is "<<real(innerC(psi,H,psi))<<endl;
    for(int n=1;n<=nsw;n++)
    {
        if(n<n_GSE)
        {
            vector<Real> epsilonK;
            for(int i=1;i<=Korder-1;i++)
            {
                epsilonK.push_back(pow(0.5,i));
            }
            addBasis(psi,H,epsilonK,{"Cutoff",1E-8,"Method","DensityMatrix","KrylovOrd",Korder,"DoNormalize",true,"NumCenter",1});
        }
        auto energy=tdvp(psi,H,-t,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"ErrGoal",1E-8,"NumCenter",2});
        psi.normalize();
        cout<<t*n<<" "<<real(innerC(psi,H,psi))<<endl;
    }
    ofstream fp;
    fp.open("./data/CorrelatorsN"+to_string(N)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(h)+".txt");
    for(int k=1;k<=N-1;k++)
    {
        auto XX=correlators2point(0,k,"Sx","Sx",psi);

    }
    fp.close();
}
int main(int argc, char* argv[])
{
    Args::global().add("Cutoff", 1E-10);
    Args::global().add("Method","DensityMatrix");
    Args::global().add("KrylovOrd",5);
    Args::global().add("DoNormalize",true);
    Args::global().add("MaxDim",5000);
    double alpha=stod(argv[1]);
    N=stoi(argv[2]);
    sites=SpinHalf(N,{"ConserveQNs=",false});
    double t=stod(argv[3]);
    double tend=stod(argv[4]);
    int n_GSE=stoi(argv[5]);
    int Korder=stoi(argv[6]);
    double h=stod(argv[7]);
    imaginarytimeevolutionTDVP(alpha,t,tend,n_GSE,Korder,h);
}