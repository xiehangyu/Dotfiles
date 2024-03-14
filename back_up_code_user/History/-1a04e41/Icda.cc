#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
using namespace itensor;

double Period_length(int i, int j,int N)
{//The spins are in a 1D chain, so the distance between them is either i-j or N-(i-j)
    return std::min(fabs(i-j+0.0),N-fabs(i-j+0.0));
   // return fabs(i-j+0.0);
}
MPO totalspin(SpinHalf sites, int N)
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
    auto O=toMPO(ampo,{"Exact=",true});
    printf("Maximum bond dimension of O is %d",maxLinkDim(O));
    return O;
}
int 
main(int argc, char* argv[])
    {
	Real alpha = atof(argv[1]);
	auto t = atof(argv[2]);
	Real tend = atof(argv[3]);
	int nsw = tend/t;
	int dk = 3;

	printfln("alpha=%.5f",alpha);
	printfln("tstep=%.5f",t);
	printfln("tend=%.5f",tend);

	int N = atoi(argv[4]);
	printfln("N=%d",N);

    int ngse = atoi(argv[5]);

	auto sites = SpinHalf(N);
	println("Writing sites to disk");
    std::string sitesdir = "sites/";
	writeToFile(sitesdir+"alpha"+std::string(argv[1])+"size"+std::string(argv[4])+"sites",sites);

	auto ampo = AutoMPO(sites);
	double normfactor=1.0/pow(N,1-alpha);
    // Notice here each term will be counted twice
	for(int i = 1; i <= N; ++i)
        for(int j = 1+i; j <= N; ++j)
            {
			Real decayij = normfactor/std::pow(Period_length(i,j,N)+1.0,alpha);
			ampo += decayij,"Sz",i,"Sz",j;
		    }
    
	auto H = toMPO(ampo);
	printfln("Maximum bond dimension of H is %d",maxLinkDim(H));

	auto is = inds(sites);
	for(int i = 1; i <= length(is); ++i)
		is(i) = removeQNs(is(i));
	auto sitesQN = SpinHalf(is);

	ampo = AutoMPO(sitesQN);
	for(int i = 1; i <= N; ++i)
		{
		for(int j = 1; j <= N; ++j)
			{
			ampo +=  0.25,"S+",i,"S+",j;
			ampo +=  0.25,"S+",i,"S-",j;
			ampo +=  0.25,"S-",i,"S+",j;
			ampo +=  0.25,"S-",i,"S-",j;
			}
		}
	auto SxSx = toMPO(ampo);
	auto spintotal=totalspin(sites,N)
 
	ampo = AutoMPO(sitesQN);
	for(int i = 1; i <= N; ++i)
		{
		for(int j = 1; j <= N; ++j)
			{
			ampo +=  -0.25,"S+",i,"S+",j;
			ampo +=   0.25,"S+",i,"S-",j;
			ampo +=   0.25,"S-",i,"S+",j;
			ampo +=  -0.25,"S-",i,"S-",j;
			}
		}
	auto SySy = toMPO(ampo);

	ampo = AutoMPO(sitesQN);
	for(int i = 1; i <= N; ++i)
		{
		for(int j = 1; j <= N; ++j)
			{
			ampo += "Sz",i,"Sz",j;
			}
		}
	auto SzSz = toMPO(ampo);

	ampo = AutoMPO(sitesQN);
	for(int i = 1; i <= N; ++i)
		{
		for(int j = 1; j <= N; ++j)
			{
			ampo +=  -0.5*1_i,"S+",i,"Sz",j;
			ampo +=   0.5*1_i,"S-",i,"Sz",j;
			}
		}
	auto SySz = toMPO(ampo);

	ampo = AutoMPO(sitesQN);
	for(int i = 1; i <= N; ++i)
		{
		for(int j = 1; j <= N; ++j)
			{
			ampo +=  -0.5*1_i,"Sz",i,"S+",j;
			ampo +=   0.5*1_i,"Sz",i,"S-",j;
			}
		}
	auto SzSy = toMPO(ampo);

	auto state = InitState(sites);
	for(int i = 1; i <= N; ++i) 
        {
		state.set(i,"Up");
		}
	auto psi1 = MPS(state);

	for(int i = 1; i <= N; ++i)
        {
		psi1.position(i);

		auto ind = removeQNs(sites(i));
		auto indP = prime(ind);
		auto R1 = ITensor(ind, indP);
		R1.set(ind(1),indP(1), 1./sqrt(2));
		R1.set(ind(2),indP(1),-1./sqrt(2));
		R1.set(ind(1),indP(2), 1./sqrt(2));
		R1.set(ind(2),indP(2), 1./sqrt(2));
		
		psi1.Aref(i) = noPrime(R1*psi1.A(i));
		auto si = removeQNs(sites(i));
		auto Sxi = op(SpinHalfSite(si),"Sx");
		auto Syi = op(SpinHalfSite(si),"Sy");
		Real sxi = std::real((psi1.A(i) * Sxi * dag(prime(psi1.A(i),"Site"))).cplx());
		Real syi = std::real((psi1.A(i) * Syi * dag(prime(psi1.A(i),"Site"))).cplx());
		Real szi = std::real((psi1.A(i) * sites.op("Sz",i) * dag(prime(psi1.A(i),"Site"))).cplx());
		printfln("Initial sx at site %d =  %.10f,    sy at site %d =  %.10f,    sz at site %d = %.10f", i,sxi,i,syi,i,szi);
        }
	psi1.position(1);
	printfln("norm = %.5f", real(innerC(psi1,psi1)));

	auto energy = real(innerC(psi1,H,psi1));
	printfln("Initial energy = %.5f", energy);
	
	println("----------------------------------------GSETDVP1---------------------------------------");
	
	auto sweeps = Sweeps(1);
	sweeps.maxdim() = 10000;
	sweeps.cutoff() = 1E-10;
	sweeps.niter() = 50;

	auto args = Args("ErrGoal=",1E-8,"IsHermitian=",true,"Quiet=",true,"NumCenter=",1,"Truncate",true);
    std::string tempdir = std::string("temp/")+std::string("alpha")+argv[1]+std::string(argv[4]);
	auto args_addbasis = Args("Cutoff=",1E-4,"KrylovOrd",dk,"DoNormalize",true,"NumCenter=",1,"WriteDir",tempdir,"WriteDim",1000);
 	
	auto sweeps_plain = Sweeps(1);
	sweeps_plain.maxdim() = 10000;
	sweeps_plain.cutoff() = 1E-8;
	sweeps_plain.niter() = 50;
	
	auto args_plain = Args("ErrGoal=",1E-7,"IsHermitian=",true,"Quiet=",true,"NumCenter=",2,"WriteDir",tempdir,"WriteDim",5000);

	std::string datadir = "data/";
	std::ofstream myfilex;
	std::ofstream myfileyy;
	std::ofstream myfilezz;
	std::ofstream myfileyz;
	std::ofstream myfilezy;
	std::ofstream myfilesq;
	std::ofstream myfiles2;
	myfilex.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"sx"+std::string(argv[5])+"x"+std::string(argv[5]));
	myfileyy.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"sysy"+std::string(argv[5])+"x"+std::string(argv[5]));
	myfilezz.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"szsz"+std::string(argv[5])+"x"+std::string(argv[5]));
	myfileyz.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"sysz"+std::string(argv[5])+"x"+std::string(argv[5]));
	myfilezy.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"szsy"+std::string(argv[5])+"x"+std::string(argv[5]));
	myfilesq.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"squeez"+std::string(argv[5])+"x"+std::string(argv[5]));
	myfiles2.open(datadir+std::string("alpha")+argv[1]+std::string("delta")+std::string(argv[2])+"s2"+std::string(argv[5])+"x"+std::string(argv[5]));

	auto sx = 0.0;
	for(int j = 1; j <= N; ++j)
        	{
        	psi1.position(j);
		auto Sxj = op(SpinHalfSite(removeQNs(sites(j))),"Sx");
		sx += std::real((psi1.A(j) * Sxj * dag(prime(psi1.A(j),"Site"))).cplx());
		}
	psi1.position(1);
	auto sxsx = real(innerC(psi1,SxSx,psi1));
	auto sysy = real(innerC(psi1,SySy,psi1));
	auto szsz = real(innerC(psi1,SzSz,psi1));
	auto sysz = real(innerC(psi1,SySz,psi1));
	auto szsy = real(innerC(psi1,SzSy,psi1));
	auto squeez = -10*std::log10(std::pow(std::sqrt(N)*std::sqrt(0.5*(sysy+szsz)-0.5*std::sqrt(std::pow(sysy-szsz,2)+std::pow(sysz+szsy,2)))/std::abs(sx),2));

	myfilex << std::setprecision(10)<< sx<< std::endl;
	myfileyy << std::setprecision(10)<< sysy<< std::endl;
	myfilezz << std::setprecision(10)<< szsz<< std::endl;
	myfileyz << std::setprecision(10)<< sysz<< std::endl;
	myfilezy << std::setprecision(10)<< szsy<< std::endl;
	myfilesq << std::setprecision(10)<< squeez<< std::endl;
	myfiles2 << std::setprecision(10)<< sxsx+sysy+szsz<< std::endl;

	std::string psidir = "psi/";
	LocalMPO PH(H,args_plain);
	for(int sw = 1; sw <= nsw; ++sw)
		{
		if(sw < ngse)
		    {
		    std::vector<Real> trunc = {1E-8,1E-6};
		    addBasis(psi1, H, trunc, args_addbasis);
		    energy = tdvp(psi1,H,-1_i*t,sweeps,args);
		    }
		else
		    {
	        printfln("Read MPS from sweep %d", sw-1); 
            psi1 = readFromFile<MPS>(psidir+"alpha"+std::string(argv[1])+"psi"+std::string(argv[4])+"sweep"+std::to_string(sw-1));
		    energy = TDVPWorker(psi1,PH,-1_i*t,sweeps_plain,args_plain);
		    //energy = tdvp(psi1,H,-1_i*t,sweeps_plain,args_plain);
		    }
		println("Writing wavefunction 'psi' to disk");
		writeToFile(psidir+"alpha"+std::string(argv[1])+"psi"+std::string(argv[4])+"sweep"+std::to_string(sw),psi1);

		sx = 0.0;
		for(int j = 1; j <= N; ++j)
			{
			psi1.position(j);
			auto Sxj = op(SpinHalfSite(removeQNs(sites(j))),"Sx");
			sx += std::real((psi1.A(j) * Sxj * dag(prime(psi1.A(j),"Site"))).cplx());
			}
		psi1.position(1);//position does not preserve linkindex of MPS
		sxsx = real(innerC(psi1,SxSx,psi1));
		sysy = real(innerC(psi1,SySy,psi1));
		szsz = real(innerC(psi1,SzSz,psi1));
		sysz = real(innerC(psi1,SySz,psi1));
		szsy = real(innerC(psi1,SzSy,psi1));
		squeez = -10*std::log10(std::pow(std::sqrt(N)*std::sqrt(0.5*(sysy+szsz)-0.5*std::sqrt(std::pow(sysy-szsz,2)+std::pow(sysz+szsy,2)))/std::abs(sx),2));
	
		myfilex << std::setprecision(10)<< sx<< std::endl;
		myfileyy << std::setprecision(10)<< sysy<< std::endl;
		myfilezz << std::setprecision(10)<< szsz<< std::endl;
		myfileyz << std::setprecision(10)<< sysz<< std::endl;
		myfilezy << std::setprecision(10)<< szsy<< std::endl;
		myfilesq << std::setprecision(10)<< squeez<< std::endl;
		myfiles2 << std::setprecision(10)<< sxsx+sysy+szsz<< std::endl;
        printfln("Measurement at sweep %d finished", sw);
		std::cout<<"The expectation of spins are "<<N/2*(N/2+1)-real(innerC(psi1,spintotal,psi1))<<endl;
		}
	myfilex.close();
	myfileyy.close();
	myfileyz.close();
	myfilezy.close();
	myfilesq.close();
	myfiles2.close();

	return 0;
    }
