#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace itensor;
using namespace std;

//N is the number of sites, input by the user
static int N;
static SiteSet sites;

double Period_length(int i, int j)
{//The spins are in a 1D chain, so the distance between them is either i-j or N-(i-j)
//    return min(fabs(i-j+0.0),N-fabs(i-j+0.0));
    return fabs(i-j+0.0);
}

auto LongrangeHamiltonian(double alpha,double J)
{//Construct the long-range Hamiltonian with exponent alpha
    auto ampo=AutoMPO(sites);
    double normfactor=1.0/pow(N,1-alpha);
    for(int i=1;i<=N;i++)
    {
        for(int j=i+1;j<=N;j++)
        {
            ampo+=4.0*normfactor/pow(1.0+Period_length(i,j),alpha),"Sz",i,"Sz",j;
        }
        ampo += 2.0*J,"Sx",i;
    }
    auto H=toMPO(ampo);
    printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    return H;
}
auto LongrangeHamiltonianampo(double alpha,double J)
{//Construct the long-range Hamiltonian with exponent alpha
    auto ampo=AutoMPO(sites);
    double normfactor=1.0/pow(N,1-alpha);
    for(int i=1;i<=N;i++)
    {
        for(int j=i+1;j<=N;j++)
        {
            ampo+=4.0*normfactor/pow(1.0+Period_length(i,j),alpha),"Sz",i,"Sz",j;
        }
        ampo += 2.0*J,"Sx",i;
    }
    return ampo;
}
auto NearestHamiltonian(double J)
{//Construct the long-range Hamiltonian with exponent alpha
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=N-1;i++)
    {
        ampo+=4.0,"Sz",i,"Sz",i+1;
        ampo += 2.0*J,"Sx",i;
    }
    auto H=toMPO(ampo);
    printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    return H;
}

auto randomHamiltonian(double J)
{
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=N;i++)
    {
        ampo+=J*(((double)rand())/(RAND_MAX+0.0)-0.5),"Sz",i;
    }
    auto H=toMPO(ampo);
    printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    return H;
}


double vn_entropy(MPS psi, int site)
{//Calculate the von Neumann entropy of the subsystem from site 1 to site "site", the log bases is 2
    psi.position(site);
    auto l=leftLinkIndex(psi,site);
    auto s=siteIndex(psi,site);
    auto [U,S,V]=svd(psi(site),{l,s});
    auto u=commonIndex(U,S);
    Real SvN=0;
    for(auto p : range1(dim(u)))
    {
        auto lambda=elt(S,p,p);
        auto sn = sqr(lambda);
        if(sn>1E-12) SvN+=-sn*log2(sn);
    }
    return SvN;
}

auto productstate()
{//construct the initial state
    auto state=InitState(sites,"Up");
    return MPS(state);
}


MPO totalX()
{//construct the total spin X, this function is not used in the program
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=N;i++)
    {
        ampo+=1.0,"Sx",i;
    }
    return toMPO(ampo);
}

MPO totalY()
{//construct the total spin Y, this function is not used in the program
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=N;i++)
    {
        ampo+=1.0,"Sy",i;
    }
    return toMPO(ampo);
}

MPO totalZ()
{//construct the total spin Z, this function is not used in the program
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=N;i++)
    {
        ampo+=1.0,"Sz",i;
    }
    return toMPO(ampo);
}

MPO totalspin()
{//construct the total spin operator
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=N;i++)
    {
        for(int j=i+1;j<=N;j++)
        {
            ampo+=1.0,"S+",i,"S-",j;
            ampo+=1.0,"S-",i,"S+",j;
            ampo+=2.0,"Sz",i,"Sz",j;
        }
        ampo += 0.75,"Id",i;
    }
    auto O=toMPO(ampo);
    printf("Maximum bond dimension of O is %d",maxLinkDim(O));
    return O;
}

void direct_application(double alpha, double t=0.1, double tend=0.5, int n_GSE=3, int Korder=3, double randomin=0,double J=0)
{
    ofstream fp;
    fp.open("./data/DirectOBCEntropyresultN"+to_string(N)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(J)+".txt");
    auto psi1=productstate();
    auto ampo=AutoMPO(sites);
    if(alpha>-0.1)
	{
        ampo=LongrangeHamiltonianampo(alpha,J);
	}
    else
    {
        cout<<"We haven't done this yet"<<endl;
        exit(1);
    }
    auto observable=totalspin();
    int nsw=round(tend/t);
    auto expH=toExpH(ampo,1_i*t);
    auto args=Args("Method=","DensityMatrix","Cutoff=",1E-12,"MaxDim=",10000);
    cout<<"The initial deviated expectation value is "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
    cout<<"The half chain entropy is "<<vn_entropy(psi1,N/2)<<endl;
    for(int n = 1; n <= nsw; ++n)
        {
            psi1=applyMPO(expH,psi1,args);
            psi1.noPrime().normalize();
        cout<<t*n<<" "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
        cout<<t*n<<" "<<vn_entropy(psi1,N/2)<<endl;
        printfln("\nMaximum bond dimension at time %.1f is %d",t*n,maxLinkDim(psi1));
        fp<<t*n<<" "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<"  "<<vn_entropy(psi1,N/2)<<endl;
	fp.flush();
        }
}
void overlap(double alpha, double t=0.1, double tend=0.5, int n_GSE=3, int Korder=3, double randomin=0,double J=0)
{//Time evolution using TDVP, the decay exponent of the interaction is alpha, the time step is t, the total time is tend, the number of global subspace expansion is n_GSE, the order of Krylov subspace is Korder, the two printed results are the deviation of the expectation value of the total spin from its maximum value and the half chain entropy.
    ofstream fp,gp;
    fp.open("./data/XXOBCEntropyresultN"+to_string(N)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(J)+".txt");
    gp.open("./data/XXOBCOverlapN"+to_string(N)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(J)+".csv");
    auto psi1=productstate();
    auto vector<MPS> psis;
    auto temp=psi1;
    psis.push_back(temp);
    temp_overlap=innerC(psi1,psi1);
    gp<<real(temp_overlap)<<","<<imag(temp_overlap)<<endl;
    auto H=MPO();
    if(alpha>-0.1)
	{
		H=LongrangeHamiltonian(alpha,J);
	}
    else
    {
	H=NearestHamiltonian(J);
    }
    auto observable=totalspin();
    int nsw=round(tend/t);
    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 10000;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() =50;
    cout<<endl;
    cout<<"The initial deviated expectation value is "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
    cout<<"The half chain entropy is "<<vn_entropy(psi1,N/2)<<endl;
    for(int n = 1; n <= nsw; ++n)
        {
        if(n < n_GSE)
            {
            // Global subspace expansion
            vector<Real> epsilonK;
            for(int i=1;i<=Korder-1;i++)
            {
                epsilonK.push_back(1E-12);
            }
            
            addBasis(psi1,H,epsilonK,{"Cutoff",1E-6,"Method","DensityMatrix","KrylovOrd",Korder,"DoNormalize",true,"NumCenter",1});
 //addBasis(psi1,H,epsilonK,{"Cutoff",1E-6,"Method","DensityMatrix","KrylovOrd",Korder,"DoNormalize",true});
            }
        
        // TDVP sweep
        
        tdvp(psi1,H,-t*1_i,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"ErrGoal",1E-8,"IsHermitian",true,"NumCenter",2});
        temp=psi1;
        psis.push_back(temp);
        
        //tdvp(psi1,H,-t*1_i,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"NumCente", 1,"ErrGoal",1E-9,"IsHermitian",true});
        if(n%10==0 && fabs(randomin)>=1E-8)
        {
           ; 
            //auto randomH=randomHamiltonian(randomin);
            //tdvp(psi1,randomH,-t*1_i,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"NumCente", 1,"ErrGoal",1E-9,"IsHermitian",true});
        }

        cout<<t*n<<" "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
        cout<<t*n<<" "<<vn_entropy(psi1,N/2)<<endl;
        fp<<t*n<<" "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<"  "<<vn_entropy(psi1,N/2)<<endl;
	    fp.flush();
        int temp_size=psis.size();
        for(state_index=0;state_index<temp_size;state_index++)
        {
            auto temp_overlap=innerC(psis[state_index],psis[temp_size-1]);
            gp<<real(temp_overlap)<<","<<imag(temp_overlap);
            if(state_index<temp_size-1)
            {
                gp<<",";
            }
        }
        gp<<endl;
        gp.flush();
        }
    fp.close();
    gp.close();
   
}
void timeevolutionTDVP(double alpha, double t=0.1, double tend=0.5, int n_GSE=3, int Korder=3, double randomin=0,double J=0)
{//Time evolution using TDVP, the decay exponent of the interaction is alpha, the time step is t, the total time is tend, the number of global subspace expansion is n_GSE, the order of Krylov subspace is Korder, the two printed results are the deviation of the expectation value of the total spin from its maximum value and the half chain entropy.
    ofstream fp;
    fp.open("./data/XXOBCEntropyresultN"+to_string(N)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(J)+".txt");
    auto psi1=productstate();
    auto H=MPO();
    if(alpha>-0.1)
	{
		H=LongrangeHamiltonian(alpha,J);
	}
    else
    {
	H=NearestHamiltonian(J);
    }
    auto observable=totalspin();
    int nsw=round(tend/t);
    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 10000;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() =50;
    cout<<endl;
    cout<<"The initial deviated expectation value is "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
    cout<<"The half chain entropy is "<<vn_entropy(psi1,N/2)<<endl;
    for(int n = 1; n <= nsw; ++n)
        {
        if(n < n_GSE)
            {
            // Global subspace expansion
            vector<Real> epsilonK;
            for(int i=1;i<=Korder-1;i++)
            {
                epsilonK.push_back(1E-12);
            }
            
addBasis(psi1,H,epsilonK,{"Cutoff",1E-6,"Method","DensityMatrix","KrylovOrd",Korder,"DoNormalize",true,"NumCenter",1});
 //addBasis(psi1,H,epsilonK,{"Cutoff",1E-6,"Method","DensityMatrix","KrylovOrd",Korder,"DoNormalize",true});
            }
        
        // TDVP sweep
        
        tdvp(psi1,H,-t*1_i,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"ErrGoal",1E-8,"IsHermitian",true,"NumCenter",2});
        
        //tdvp(psi1,H,-t*1_i,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"NumCente", 1,"ErrGoal",1E-9,"IsHermitian",true});
        if(n%10==0 && fabs(randomin)>=1E-8)
        {
           ; 
            //auto randomH=randomHamiltonian(randomin);
            //tdvp(psi1,randomH,-t*1_i,sweeps,{"Truncate",true,"DoNormalize",true,"Quiet",true,"NumCente", 1,"ErrGoal",1E-9,"IsHermitian",true});
        }

        cout<<t*n<<" "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
        cout<<t*n<<" "<<vn_entropy(psi1,N/2)<<endl;
        fp<<t*n<<" "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<"  "<<vn_entropy(psi1,N/2)<<endl;
	fp.flush();
        }
   
}

int main(int argc, char* argv[])
{//The input parameters are alpha, number of sites, time step, total time, number of global subspace expansion, order of Krylov subspace
    Args::global().add("Cutoff",1E-12);
    Args::global().add("Method","DensityMatrix");
    Args::global().add("KrylovOrd",5);
    Args::global().add("DoNormalize",true);
    Args::global().add("MaxDim",5000);
    double alpha=stod(argv[1]);
    N=stoi(argv[2]);
    sites=SpinHalf(N,{"ConserveQNs=",false});
    if(argc==3)
    {//Use the default value of the function timeevolutionTDVP
        timeevolutionTDVP(alpha);
    }
    else
    {//Use the input value from the user
        double t=stod(argv[3]);
        double tend=stod(argv[4]);
        int GSE=stoi(argv[5]);
        int Korder=stoi(argv[6]);
        double J=stod(argv[7]);
	double h=stod(argv[8]);
        direct_application(alpha,t,tend,GSE,Korder,J,h);
    }
    return 0;
}
