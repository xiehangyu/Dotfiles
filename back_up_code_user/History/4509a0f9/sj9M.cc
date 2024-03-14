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

//L is the number of sites along one dimension input by the user
static int L;
static SiteSet sites;


inline double D2distance(int i, int j)
{
    int x1=i%L;
    int y1=i/L;
    int x2=j%L;
    int y2=j/L;
    return pow((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2),0.5);
}

double Period_length(int i, int j)
{//The spins are in a 1D chain, so the distance between them is either i-j or N-(i-j)
//    return min(fabs(i-j+0.0),N-fabs(i-j+0.0));
    return fabs(i-j+0.0);
}

auto LongrangeHamiltonian(double alpha,double J)
{//Construct the long-range Hamiltonian with exponent alpha
    auto ampo=AutoMPO(sites);
    double normfactor=1.0/pow(L,2-alpha);
    for(int i=1;i<=L*L;i++)
    {
        for(int j=i+1;j<=L*L;j++)
        {
            ampo+=4*normfactor/pow(1.0+D2distance(i,j),alpha),"Sx",i,"Sx",j;
        }
        ampo += J*2,"Sz",i;
    }
    auto H=toMPO(ampo);
    printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    return H;
}

auto productstate()
{//construct the initial state
    auto state=InitState(sites,"Up");
    return MPS(state);
}


MPO totalX()
{//construct the total spin X, this function is not used in the program
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=L*L;i++)
    {
        ampo+=1.0,"Sx",i;
    }
    return toMPO(ampo);
}

MPO totalY()
{//construct the total spin Y, this function is not used in the program
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=L*L;i++)
    {
        ampo+=1.0,"Sy",i;
    }
    return toMPO(ampo);
}

MPO totalZ()
{//construct the total spin Z, this function is not used in the program
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=L*L;i++)
    {
        ampo+=1.0,"Sz",i;
    }
    return toMPO(ampo);
}

MPO totalspin()
{//construct the total spin operator
    auto ampo=AutoMPO(sites);
    for(int i=1;i<=L*L;i++)
    {
        for(int j=i+1;j<=L*L;j++)
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


void overlap(double alpha, double t=0.1, double tend=0.5, int n_GSE=3, int Korder=3, double randomin=0,double J=0)
{//Time evolution using TDVP, the decay exponent of the interaction is alpha, the time step is t, the total time is tend, the number of global subspace expansion is n_GSE, the order of Krylov subspace is Korder, the two printed results are the deviation of the expectation value of the total spin from its maximum value and the half chain entropy.
    ofstream gp;
    vector<MPS> psis;
    gp.open("./data/OBCOverlapL"+to_string(L)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(J)+".csv");
    auto psi1=productstate();
    auto temp=psi1;
    psis.push_back(temp);
    auto temp_value=innerC(psi1,psi1);
    gp<<real(temp_value)<<","<<imag(temp_value)<<endl;
    gp.flush();
    auto H=MPO();
	H=LongrangeHamiltonian(alpha,J);
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
        for(int state_index=0;state_index<temp_size;state_index++)
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
    fp.open("./data/D2OBCN"+to_string(N)+"alpha"+to_string(alpha)+"t"+to_string(t)+"tend"+to_string(tend)+"GSE"+to_string(n_GSE)+"Korder"+to_string(Korder)+"randomin"+to_string(randomin)+"J"+to_string(J)+".txt");
    auto psi1=productstate();
    auto H=MPO();
	H=LongrangeHamiltonian(alpha,J);
    auto observable=totalspin();
    int nsw=round(tend/t);
    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 10000;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() =50;
    cout<<endl;
    cout<<"The initial deviated expectation value is "<<N/2*(N/2+1)-real(innerC(psi1,observable,psi1))<<endl;
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
    L=stoi(argv[2]);
    sites=SpinHalf(L*L,{"ConserveQNs=",false});
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
        timeevolutionTDVP(alpha,t,tend,GSE,Korder,J,h);
    }
    return 0;
}
