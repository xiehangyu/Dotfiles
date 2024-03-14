#include <ctime>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/src/Core/products/Parallelizer.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <physics.hpp>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <string>
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using namespace std;
class mypower
{
    private:
        int n;
    public:
        mypower(int n=0): n(n){};
        double operator()(double x)
        {
            return pow(x,n);
        }
};        
MatrixXcd fourier_transformation_from_kspace_to_realspace_period1(int Nsites)
{
    MatrixXcd UF(Nsites,Nsites);
    for(int i=0;i<=Nsites-1;i++)
    {
        complex<double> k=2*M_PI*i/(Nsites+0.0);
        for(int j=0;j<=Nsites-1;j++)
        {
            UF(i,j)=exp(-I*k*(j+1.0));
        }
    }
    return 1/sqrt(Nsites)*UF.adjoint();
}
double entropy_calculation(double x)
{
    if (x<=0 || x>=1)
    {
        return 0;
    }
    return -x*log2(x)-(1-x)*log2(1-x);
}
double calculate_entropy(MatrixXcd covariance_matrix)
{
    covariance_matrix=(covariance_matrix+covariance_matrix.adjoint())/2;
    Eigen::SelfAdjointEigenSolver<MatrixXcd> eigenval(covariance_matrix);

    return eigenval.eigenvalues().unaryExpr([](double x){if ((x<=0)||(x>=1)) return 0.0; else return (-x*log2(x)-(1-x)*log2(1-x));}).sum();
}

double calculate_N_th_order(const int N,MatrixXcd covariance_matrix)
{
    int cols=covariance_matrix.cols();
    int rows=covariance_matrix.rows();
    MatrixXcd temp=2*covariance_matrix-MatrixXcd::Identity(rows,cols);
    temp=(temp+temp.adjoint())/2;
    Eigen::SelfAdjointEigenSolver<MatrixXcd> eigenval(temp);
    VectorXcd eigenvalues_array=eigenval.eigenvalues();
    eigenvalues_array.unaryExpr([](double x){return pow(x,4);});
    double sum=0;
    int final_index=eigenvalues_array.size();
    for (int i=0;i<=final_index-1;i++)
    {
        sum += fabs(pow(eigenvalues_array(i),2*N));
    }
    return sum;
    //return abs(MatrixXcd(temp.pow(2*N)).trace());
}

MatrixXcd covariance_matrix_at_time_t_in_real_space(double t, int Nsize)
{
    static VectorXcd energies(Nsize);
    static MatrixXcd UF(Nsize,Nsize);
    static MatrixXcd covariance_matrix_in_k_space(Nsize,Nsize);

    static int label=0;
    if (label==0)
    {
        UF=fourier_transformation_from_kspace_to_realspace_period1(Nsize);
        for(int i=0;i<=Nsize-1;i++)
        {
            double k=2*M_PI*i/(Nsize+0.0);
            energies(i)=2*cos(k);
        }
        covariance_matrix_in_k_space=MatrixXcd::Zero(Nsize,Nsize);
        for (int i=0;i<=Nsize-1;i+=2)
        {
            covariance_matrix_in_k_space(i,i)=1;
        }
        covariance_matrix_in_k_space=UF.adjoint()*covariance_matrix_in_k_space*UF;
        label=1;
    }
    MatrixXcd unitary_transofrmation;
    unitary_transofrmation=(energies*I*t).unaryExpr([](complex<double> temp){return exp(temp);}).asDiagonal();
    unitary_transofrmation=UF*unitary_transofrmation;
    return unitary_transofrmation*covariance_matrix_in_k_space*unitary_transofrmation.adjoint();
}

void Ordern_SizeN_average(int Nsize, vector<int>norders, int Samplingsize,double Initialtime, double timestep)
{
    Eigen::initParallel();
    const int order_num=norders.size();
    ofstream *fps;
    fps=new ofstream[order_num];
    double ** recordings=new double*[order_num];
    for(int i=0;i<=order_num-1;i++)
    {
        recordings[i]=new double[Nsize];
        memset(recordings[i],0,sizeof(recordings[i][0])*Nsize);
    }
    covariance_matrix_at_time_t_in_real_space(Initialtime, Nsize);
    #pragma omp parallel shared(recordings) num_threads(4)
    {
        Eigen::setNbThreads(2);
        double **self_recordings=new double*[order_num];
        for(int i=0;i<=order_num-1;i++)
        {
            self_recordings[i]=new double[Nsize];
            memset(self_recordings[i],0,sizeof(self_recordings[i][0])*Nsize);
        }
        #pragma omp for
        for(int i=1;i<=Samplingsize;i++)
        {
            if(i%100==0)
            {
                cout<<i<<endl;
            }
            double time=Initialtime+i*timestep;
            MatrixXcd real_CM(Nsize,Nsize);
            real_CM=covariance_matrix_at_time_t_in_real_space(time,Nsize);
            for(int j=1;j<=Nsize;j++)
            {
                MatrixXcd temp_matrix=real_CM.block(0,0,j,j);
                for(int k=0;k<=order_num-1;k++)
                {
                    if((norders[k]==0) && (j>=Nsize/2+1))
                    {
                        continue;
                    }
                    if( (norders[k]==0))
                    {
                        self_recordings[k][j-1]+= calculate_entropy(temp_matrix);
                        continue;
                    }
                    self_recordings[k][j-1]+=calculate_N_th_order(norders[k], temp_matrix); 
                }           
            }
        }
        #pragma omp critical
        {
            for(int k=0;k<=order_num-1;k++)
            {
                for(int j=1;j<=Nsize;j++)
                {
                    recordings[k][j-1]+=self_recordings[k][j-1];
                }
            }
        }
    }
    for(int i=0;i<=order_num-1;i++)
    {
        fps[i].open("./plot_data/Nsize"+to_string(Nsize)+"norder"+to_string(norders[i])+"Samplingsize"+to_string(Samplingsize)+"Initialtime"+to_string(Initialtime)+"timestep"+to_string(timestep)+".txt");
    }
    for(int k=0;k<=order_num-1;k++)
    {
        if (norders[k]==0)
        {
            for(int j=Nsize/2+1;j<=Nsize-1;j++)
            {
                recordings[k][j-1]=recordings[k][Nsize-j-1];
            }
            recordings[k][Nsize-1]=0;
        }
        for(int i=0;i<=Nsize-1;i++) 
        {
            fps[k]<<to_string(i+1)<<"\t"<<to_string(recordings[k][i]/(Samplingsize+0.0))<<endl;
        }
    }
    return;
}
int main()
{
    clock_t starttime,endtime;
    starttime=clock();
    vector<int>norders{0};
    Ordern_SizeN_average(100, norders, 4000, 100000000, 0.2);
    endtime=clock();
    cout<<"The total time is "<<(double)(endtime-starttime)/CLOCKS_PER_SEC<<"s"<<endl;
}