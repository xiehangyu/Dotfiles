#include <ctime>
#include <omp.h>
#include <map>
#include <complex>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/src/Core/products/Parallelizer.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <physics.hpp>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <string>
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::VectorXd;
using namespace std;
template <typename T>
using matrixself = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
static int D;
static int L;
static double x_floquet=0;
static double lambda_floquet=0;
static double onsite_random;
static string file_name;
static double alpha;
static int num_of_threads;
typedef struct temp_defi_struc
{
    int index;
    int *position;
    void self_init(int i)
    {
        index = i;
        for (int j = 0; j < D; j++)
        {
            position[j] = i % L;
            i = i / L;
        }
        return;
    }
    temp_defi_struc()
    {
        index = -3;
        position = new int[D];
        return;
    }
    temp_defi_struc &operator=(const temp_defi_struc &new_one)
    {
        index = new_one.index;
        for (int j = 0; j < D; j++)
        {
            position[j] = new_one.position[j];
        }
        return *this;
    }
    temp_defi_struc(temp_defi_struc &new_one)
    {
        index = new_one.index;
        position = new int[D];
        for (int j = 0; j < D; j++)
        {
            position[j] = new_one.position[j];
        }
    }

    ~temp_defi_struc();
} node;
temp_defi_struc::~temp_defi_struc()
{
    delete[] position;
    return;
}
static node *nodes;
static int Nsites;
double distance(node node1, node node2)
{
    double distance = 0;
    for (int i = 0; i < D; i++)
    {
        distance += pow(min(abs(node1.position[i] - node2.position[i]), L - abs(node1.position[i] - node2.position[i])), 2);
    }
    return sqrt(distance);
}
complex<double> Z_kt1(double k)
{
    static double t1=M_PI/4-x_floquet;
    complex<double> tempans=exp(complex<double>(0.0,1.0)*2.0*t1)*tan(k/2.0);
    tempans += exp(complex<double>(0.0,-1.0)*2.0*t1)/tan(k/2.0);
    return tempans/(2.0*sin(2.0*t1));
}
void a_k_floquet(double k, vector<complex<double>> & returnlist)
{
    complex<double> zkt1=Z_kt1(k);
    double denominator=1.0+abs(zkt1)*abs(zkt1);
    complex<double> temp=abs(zkt1)*abs(zkt1)*exp(-4.0*lambda_floquet)-exp(4.0*complex<double>(0.0,1.0)*x_floquet);
    returnlist[0]=temp/denominator;
    temp=(exp(-4.0*lambda_floquet)+exp(4.0*complex<double>(0.0,1.0)*x_floquet))*conj(zkt1);
    returnlist[1]=temp/denominator;
    temp=-zkt1*(exp(4.0*lambda_floquet)+exp(4.0*complex<double>(0.0,-1.0)*x_floquet));
    returnlist[2]=temp/denominator;
    temp=abs(zkt1)*abs(zkt1)*exp(4.0*lambda_floquet)-exp(4.0*complex<double>(0.0,-1.0)*x_floquet);
    returnlist[3]=temp/denominator;
    return;
}
complex<double> f_nk_floquet(double k, int n)
{
    static double measure=4.0*pow(cos(2.0*x_floquet),4)/pow(sinh(2.0*lambda_floquet),2)+2.0*cos(4.0*x_floquet);
    if(k==0||k==M_PI)
    {
        return 0;
    }
    vector<complex<double>> Mk(4);
    double temp_measure=tan(k/2.0)*tan(k/2.0);
    if(temp_measure+1.0/temp_measure>measure)
    {
        a_k_floquet(k,Mk);
        complex<double> ans1=(Mk[0]-Mk[3]+sqrt((Mk[0]+Mk[3])*(Mk[0]+Mk[3])-4.0))/(2.0*Mk[2]);
        complex<double> ans2=(Mk[0]-Mk[3]-sqrt((Mk[0]+Mk[3])*(Mk[0]+Mk[3])-4.0))/(2.0*Mk[2]);
        if(abs(Mk[2]*ans1+Mk[3])>abs(Mk[2]*ans2+Mk[3]))
        {
            return ans1;
        }
        else
        {
            return ans2;
        }
    }
    cout<<"Enter the critial region with k="<<k<<endl;
    double thetak;
    complex<double> zkt1=Z_kt1(k);
    double temp;
    temp=abs(zkt1)*abs(zkt1)*cosh(4.0*lambda_floquet)-cos(4.0*x_floquet);
    temp=temp/(1.0+abs(zkt1)*abs(zkt1));
    if((fabs(temp)>=1-1E-6) && (fabs(temp)<=1+1E-6))
    {
        thetak=0;
    }
    if(fabs(temp)>1+1E-6)
    {
        cout<<"error in acos"<<endl;
        exit(1);
    }
    if(fabs(temp)<1-1E-6)
    {
        thetak=acos(temp);
    }
    a_k_floquet(k,Mk);
    complex<double> returnvalue;
    if(thetak==0)
    {
        returnvalue=Mk[1]*(n+0.0);
        returnvalue=returnvalue/(1.0+(1.0-Mk[0])*(n+0.0));
    }
    else
    {
        returnvalue=Mk[1]*sin(n*thetak);
        returnvalue=returnvalue/((sin(thetak)*cos(n*thetak))+(cos(thetak)-Mk[0])*sin(n*thetak));
    }
    if(isnan(real(returnvalue))||isnan(imag(returnvalue)))
    {
        cout<<"nan error"<<endl;
        cout<<k<<" "<<n<<" "<<thetak<<" "<<Mk[0]<<" "<<Mk[1]<<" "<<Mk[2]<<" "<<Mk[3]<<endl;
        exit(1);
    }
    return returnvalue;
}
complex<double> f_nk_floquet2(double k, int n)
{
    MatrixXcd Mk(2,2);
    vector<complex<double>> a_k(4);
    if(k==0||k==M_PI)
    {
        Mk(0,0)=exp(-4.0*lambda_floquet);
        Mk(0,1)=0;
        Mk(1,0)=0;
        Mk(1,1)=exp(4.0*lambda_floquet);
    }
    else
    {
        a_k_floquet(k,a_k);
        Mk(0,0)=a_k[0];
        Mk(0,1)=a_k[1];
        Mk(1,0)=a_k[2];
        Mk(1,1)=a_k[3];
    }
    VectorXcd initstate(2);
    initstate(0)=0.0;
    initstate(1)=1.0;
    VectorXcd finalstate(2);
    finalstate = Mk.pow(n)*initstate;
    return(finalstate(0)/finalstate(1));
}
complex<double> varphi_j(int j, const vector<complex<double>> &fs)
{
    complex<double> sum=0;
    for(int i=0;i<=L-1;i++)
    {
        double k=2.0*M_PI*i/L;
        sum += exp(complex<double>(0.0,-1.0)*k*(j+0.0))*(fs[i]+conj(fs[i]))/(1.0+abs(fs[i])*abs(fs[i]));
    }
    sum=sum/(L+0.0);
    sum=sum*complex<double>(0.0,1.0);
    if (fabs(sum.imag())>1E-6)
    {
        cout<<"error in varphi_j"<<endl;
        exit(1);
    }
    return sum;
}
complex<double> psi_j(int j, const vector<complex<double>> &fs)
{
    complex<double> sum=0;
    for(int i=0;i<=L-1;i++)
    {
        double k=2.0*M_PI*i/L;
        double norm=abs(fs[i])*abs(fs[i]);
        sum += exp(complex<double>(0.0,-1.0)*k*(j+0.0))*(fs[i]-conj(fs[i])+norm-1.0)/(1.0+norm);
    }
    sum=sum/(L+0.0);
    if (fabs(sum.imag())>1E-6)
    {
        cout<<"error in psi_j"<<endl;
        exit(1);
    }
    return sum;
}

void generating_majoranamatrix(MatrixXcd &Majorana_matrix, int n)
{
    vector<complex<double>> fs(L);
    for(int i=0;i<L;i++)
    {
        fs[i]=f_nk_floquet(2.0*M_PI*i/L,n);
    }
    map<int,complex<double>> varphi_j_map;
    map<int,complex<double>> psi_j_map;
    #pragma omp parallel for shared(varphi_j_map,psi_j_map) num_threads(num_of_threads)
    for(int i=-L+1;i<=L-1;i++)
    {
        complex<double> temp=varphi_j(i,fs);
        complex<double> temp2=psi_j(i,fs);
        #pragma omp critical
        {
            varphi_j_map[i]=temp;
            psi_j_map[i]=temp2;
        }
    }
    for(int i=0;i<=L-1;i++)
    {
        for(int j=0;j<=L-1;j++)
        {
            Majorana_matrix(2*i,2*j)=-varphi_j_map[j-i];
            Majorana_matrix(2*i,2*j+1)=psi_j_map[j-i];
            Majorana_matrix(2*i+1,2*j)=-psi_j_map[i-j];
            Majorana_matrix(2*i+1,2*i+1)=varphi_j_map[j-i];
        }
    }
}
void verify_covariance_scaling(int n)
{

    ofstream fp;
    string filename;
    filename="./plot_data/covariance_scaling n="+to_string(n)+".txt";
    fp.open(filename);
    int old_L=L;
    if(L>20000)
    {
        L=L/2;
    }
    else
    {
        L=10000;
    }
    vector<complex<double>> fs(L);
    for(int i=0;i<=L-1;i++)
    {
        fs[i]=f_nk_floquet(2.0*M_PI*i/L,n);
    }
    for(int i=0;i<=L-1;i++)
    {
        fp<<i<<"\t"<<abs(varphi_j(i,fs))<<"\t"<<abs(psi_j(i,fs))<<"\t"<<abs(fs[i])<<"\t"<<real(fs[i])<<"\t"<<imag(fs[i])<<endl;
        //cout<<i<<"\t"<<abs(varphi_j(i,fs))<<"\t"<<abs(psi_j(i,fs))<<"\t"<<abs(fs[i])<<"\t"<<real(fs[i])<<"\t"<<imag(fs[i])<<endl;
    }
    L=old_L;
}
void initialization_node()
{
    nodes = new node[Nsites];
    for (int i = 0; i < Nsites; i++)
    {
        nodes[i].self_init(i);
    }
}
void initialization_Hamiltonian_double(MatrixXd &Hamiltonian)
{
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = i + 1; j < Nsites; j++)
        {
            Hamiltonian(i, j) = 1.0 / pow(distance(nodes[i], nodes[j]), alpha);
            int key = 0;
            for (int m = 0; m < D; m++)
            {
                key += nodes[i].position[m] - nodes[j].position[m];
            }
            key = abs(key);
            if (key % 2 == 1)
            {
                if (key % 4 == 3)
                {
                    Hamiltonian(i, j) *= +1;
                }
            }
            if (key % 2 == 0)
            {
                Hamiltonian(i, j) *= +1;
            }
            Hamiltonian(j, i) = Hamiltonian(i, j);
        }
        Hamiltonian(i, i) = onsite_random * (rand() / (RAND_MAX + 0.0) - 0.5) * 2;
    }
}
std::complex<double> Fourier_transform_imaginary_2D(node &node1, node &node2)
{
    static std::complex<double> **lists;
    static int key = 0;
    if (key == 1)
    {
        return lists[(node2.position[0] - node1.position[0] + L) % L][(node2.position[1] - node1.position[1] + L) % L];
    }
    lists = new std::complex<double> *[L];
    for (int i = 0; i < L; i++)
    {
        lists[i] = new std::complex<double>[L];
    }
#pragma omp parallel for shared(lists) num_threads(num_of_threads)
    for (int x = 0; x < L; x++)
    {
        if (x > (L + 1) / 2)
        {
            for (int y = 0; y < L; y++)
            {
                lists[x][y] = conj(lists[L - x][L - y]);
            }
            continue;
        }
        for (int y = 0; y < L; y++)
        {
            std::complex<double> sum = 0;
            double x1 = (double)(x);
            double x2 = (double)(y);
            for (int i = 0; i <= L - 1; i++)
            {
                double k1 = 2 * M_PI * i / (L + 0.0);
                for (int k = 0; k <= L - 1; k++)
                {
                    double k2 = 2 * M_PI * k / (L + 0.0);
                    if (sin(k1) * sin(k1) + sin(k2) * sin(k2) == 0)
                    {
                        continue;
                    }
                    sum += sin(k2) / sqrt(sin(k1) * sin(k1) + sin(k2) * sin(k2)) * exp(std::complex<double>(0, 1) * (k1 * x1 + k2 * x2));
                }
            }
            lists[x][y] = sum / (L * L + 0.0);
        }
    }
    key = 1;
    return lists[(node2.position[0] - node1.position[0] + L) % L][(node2.position[1] - node1.position[1] + L) % L];
}
std::complex<double> Fourier_transform_real_2D(node &node1, node &node2)
{
    static std::complex<double> **lists;
    static int key = 0;
    if (key == 1)
    {
        return lists[(node2.position[0] - node1.position[0] + L) % L][(node2.position[1] - node1.position[1] + L) % L];
    }
    lists = new std::complex<double> *[L];
    for (int i = 0; i < L; i++)
    {
        lists[i] = new std::complex<double>[L];
    }
#pragma omp parallel for shared(lists) num_threads(num_of_threads)
    for (int x = 0; x < L; x++)
    {
        if (x > (L + 1) / 2)
        {
            for (int y = 0; y < L; y++)
            {
                lists[x][y] = conj(lists[L - x][L - y]);
            }
            continue;
        }
        for (int y = 0; y < L; y++)
        {
            std::complex<double> sum = 0;
            double x1 = (double)(x);
            double x2 = (double)(y);
            for (int i = 0; i <= L - 1; i++)
            {
                double k1 = 2 * M_PI * i / (L + 0.0);
                for (int k = 0; k <= L - 1; k++)
                {
                    double k2 = 2 * M_PI * k / (L + 0.0);
                    if (sin(k1) * sin(k1) + sin(k2) * sin(k2) == 0)
                    {
                        continue;
                    }
                    sum += sin(k1) / sqrt(sin(k1) * sin(k1) + sin(k2) * sin(k2)) * exp(std::complex<double>(0, 1) * (k1 * x1 + k2 * x2));
                }
            }
            lists[x][y] = sum / (L * L + 0.0);
        }
    }
    key = 1;
    ofstream fp;
    fp.open("trial.txt");
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            fp << lists[i][j] << "\t";
        }
        fp << endl;
    }
    return lists[(node2.position[0] - node1.position[0] + L) % L][(node2.position[1] - node1.position[1] + L) % L];
}
void initialization_Hamiltonian_double_2D(MatrixXd & Hamiltonian)
{
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = i + 1; j < Nsites; j++)
        {
            int x1=nodes[j].position[0]-nodes[i].position[0];
            int x2=nodes[j].position[1]-nodes[i].position[1];
            if(x1>L/2)
            {
                x1-=L;
            }
            if(x1<-L/2)
            {
                x1+=L;
            }
            if(x2>L/2)
            {
                x2-=L;
            }
            if(x2<-L/2)
            {
                x2+=L;
            }
            if ((abs(x1)%2==1)&&(abs(x2)%2==1))
            {
                Hamiltonian(i,j)=8.0/(x1*x2+0.0);
            }
            else
            {
                Hamiltonian(i,j)=0;
            }
            Hamiltonian(j,i)=Hamiltonian(i,j);
        }
        Hamiltonian(i, i) = onsite_random * (rand() / (RAND_MAX + 0.0) - 0.5) * 2;
    }
    //cout<<Hamiltonian<<endl;
}
void initialization_Hamiltonian_1D_SSH2(MatrixXd &Hamiltonian)
{
    /*
    a_m=\sum_k a_k e^(-ikm)/sqrt(L)
    */
   double value;
   double distance_value;
   for(int i=0;i<Nsites;i+=2)
   {
        for(int j=i+2;j<Nsites;j+=2)
        {
            distance_value=min(fabs(j-i),L-fabs(j-i))/2.0;
            value=(cos(M_PI*(distance_value)))/pow(distance_value,alpha);
            Hamiltonian(i,j)=value;
            Hamiltonian(j,i)=value;
        }
   }
   for(int i=1;i<Nsites;i+=2)
   {
        for(int j=i+2;j<Nsites;j+=2)
        {
            distance_value=min(fabs(j-i),L-fabs(j-i))/2.0;
            value=-(cos(M_PI*(distance_value)))/pow(distance_value,alpha);
            Hamiltonian(i,j)=value;
            Hamiltonian(j,i)=value;
        }
   }
   for(int i=0;i<Nsites;i+=2)
   {
        for(int j=1;j<Nsites;j+=2)
        {
            if(((j-i<=L/2)&&(j-i>0))||(i-j>=L/2))
            {
                if((j-i==1)||(j-i==-L+1))
                {
                    continue;
                }
                distance_value=min(fabs(j-i),L-fabs(j-i));
                distance_value=(distance_value-1)/2.0;
                value=1/pow(distance_value,alpha);
            }
            else
            {
                distance_value=min(fabs(j-i),L-fabs(j-i));
                distance_value=(distance_value+1)/2.0;
                value=-1/pow(distance_value,alpha);
            }
            Hamiltonian(i,j)=value;
            Hamiltonian(j,i)=value;
        }
   }
}
void initialization_Hamiltonian_1D_SSH(MatrixXd &Hamiltonian)
{
    double value;
    double distance_value;
    for(int i=0;i<Nsites;i+=2)
    {
        for(int j=1;j<Nsites;j+=2)
        {
            if(((j-i<=L/2)&&(j-i>0))||(i-j>=L/2))
            {
                distance_value=min(fabs(j-i),L-fabs(j-i));
                distance_value=(distance_value+1)/2;
                value=(cos(M_PI*(distance_value))+1)/pow(distance_value,alpha);
            }
            else
            {
                distance_value=min(fabs(j-i),L-fabs(j-i));
                distance_value=(distance_value+1)/2;
                value=(cos(M_PI*(distance_value))-1)/pow(distance_value,alpha);
            }
            Hamiltonian(i,j)=value;
            Hamiltonian(j,i)=value;
        }
    }
}
void initialization_Hamiltonian_complex(MatrixXcd &Hamiltonian)
{
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = i + 1; j < Nsites; j++)
        {
            std::complex<double> temp1;
            temp1 = Fourier_transform_real_2D(nodes[i], nodes[j]);
            std::complex<double> temp2;
            temp2 = Fourier_transform_imaginary_2D(nodes[i], nodes[j]);
            Hamiltonian(i, j) = temp1 + temp2 * std::complex<double>(0, 1);
            Hamiltonian(j, i) = temp1 - temp2 * std::complex<double>(0, 1);
        }
        Hamiltonian(i, i) = onsite_random * (rand() / (RAND_MAX + 0.0) - 0.5) * 2;
    }
}
template <typename T>
matrixself<T> ground_state(matrixself<T> &Hamiltonian)
{
    int i;
    Hamiltonian = (Hamiltonian + Hamiltonian.adjoint()) / 2;
    Eigen::SelfAdjointEigenSolver<matrixself<T>> es;
    es.compute(Hamiltonian);
    VectorXd energies = es.eigenvalues();
    matrixself<T> eigenvectors = es.eigenvectors();
    matrixself<T> groundstate = matrixself<T>::Zero(Nsites, Nsites);
    for (i = 0; i < Nsites; i++)
    {
        if (energies[i] < 0)
        {
            groundstate(i, i) = 1;
        }
        else
        {
            cout << "The gap is " << energies[i] << endl;
            break;
        }
    }
    ofstream fp;
    fp.open("./plot_data/" + file_name + "EnergyD" + to_string(D) + "L" + to_string(L) + "alpha" + to_string(alpha) + "onsiterandom" + to_string(onsite_random) + ".txt");
    fp << energies.transpose() << endl;
    fp.close();
    return eigenvectors.adjoint().transpose() * groundstate * eigenvectors.transpose();
}
MatrixXcd normalized_gamma_matrix(MatrixXcd &gamma_matrix)
{
    MatrixXcd normalized_matrix;
    normalized_matrix = complex<double>(0.0,1.0)*gamma_matrix;
    for(int i=0;i<Nsites;i++)
    {
        normalized_matrix(i,i)=0;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es;
    normalized_matrix = (normalized_matrix + normalized_matrix.adjoint()) / 2;
    es.compute(normalized_matrix);
    VectorXd energies = es.eigenvalues();
    MatrixXcd eigenvectors = es.eigenvectors();
    normalized_matrix = MatrixXcd::Zero(Nsites, Nsites);
    for (int i = 0; i < Nsites; i++)
    {
        if (energies[i] < 0)
        {
            normalized_matrix(i , i) = -1;
        }
        if (energies[i] > 0)
        {
            normalized_matrix(i , i) = 1;
        }
    }
    return eigenvectors * normalized_matrix * eigenvectors.adjoint();
}
template <typename T>
double entropy_calculation(matrixself<T> &covariance_matrix)
{
    cout << covariance_matrix.rows() << "\t" << covariance_matrix.cols() << endl;
    Eigen::SelfAdjointEigenSolver<matrixself<T>> es;
    covariance_matrix = (covariance_matrix + covariance_matrix.adjoint()) / 2;
    es.compute(covariance_matrix);
    return es.eigenvalues().unaryExpr([](double x)
                                      {if ((x<=0)||(x>=1)) return 0.0; else return(-x*log2(x)-(1-x)*log2(1-x)); })
        .sum();
}
double entropy_calculation_from_MajoranaGamma(MatrixXcd & gamma_matrix)
{
    cout<< gamma_matrix.rows() << "\t" << gamma_matrix.cols() << endl;
    MatrixXcd cova=complex<double>(0,1)*gamma_matrix;
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es;
    cova=(cova+cova.adjoint())/2;
    es.compute(cova);
    if(gamma_matrix.rows()==L)
    {
        cout<<es.eigenvalues()<<endl;
    }
    return es.eigenvalues().unaryExpr([](double x){if ((x<=-1.0)||(x>=1.0)) return 0.0; else return(-((1.0-x)/2.0)*log2((1.0-x)/2.0)); }).sum();
}
inline int whether_in_range(int start, int end, int index)
{
    if (end < L)
    {
        return (index >= start) && (index < end);
    }
    else
    {
        return (index >= start) || (index < end - L);
    }
}
template <typename T>
matrixself<T> covariance_matrix_slice(matrixself<T> &covariance_matrix, int slice, int *startpoint)
{
    vector<int> indexs;
    int *positions = new int[D];
    int total_num = 1;
    int index;
    int temp_multi;
    int temp_i;
    for (int i = 0; i < D; i++)
    {
        total_num *= slice;
    }
    for (int i = 0; i <= total_num - 1; i++)
    {

        temp_i = i;
        for (int j = 0; j < D; j++)
        {
            positions[j] = (temp_i % slice + startpoint[j]) % L;
            temp_i = temp_i / slice;
        }
        index = 0;
        temp_multi = 1;
        for (int j = 0; j < D; j++)
        {
            index += positions[j] * temp_multi;
            temp_multi *= L;
        }
        indexs.push_back(index);
    }
    delete[] positions;
    return covariance_matrix(indexs, indexs);
}
void verify_entropy_simple(int slice_total)
{
    double entropy;
    ofstream fp;
    fp.open("./plot_data/" + file_name + "EntropyD" + to_string(D) + "L" + to_string(L) + "alpha" + to_string(alpha) + "uniform" + to_string(onsite_random) + ".txt");
    MatrixXd Hamiltonian=MatrixXd::Zero(Nsites,Nsites);
    initialization_node();
    initialization_Hamiltonian_1D_SSH2(Hamiltonian);
    auto groundstate=ground_state(Hamiltonian);
    int * startpoing;
    startpoing=new int[D]{};
    #pragma omp parallel for private(entropy) schedule(dynamic) num_threads(num_of_threads/2)
    for(int i=1;i<=slice_total;i++)
    {
        auto covariance_matrix=covariance_matrix_slice(groundstate,i,startpoing);
        entropy=entropy_calculation(covariance_matrix);
        #pragma omp critical
        {
        fp<<i<<"\t"<<entropy<<endl;
        }
    }
    fp.close();    
}

void verify_entropy_simpleH(int slice_total, MatrixXcd & Hamiltonian)
{
    double entropy;
    ofstream fp;
    fp.open("./plot_data/" + file_name + "EntropyD" + to_string(D) + "L" + to_string(L) + "alpha" + to_string(alpha) + "uniform" + to_string(onsite_random) + ".txt");
    cout<<"Now calculate ground state"<<endl;
    auto groundstate=ground_state(Hamiltonian);
    cout<<"Now calculate entropy"<<endl;
    int * startpoing;
    startpoing=new int[D]{};
    #pragma omp parallel for private(entropy) schedule(dynamic) num_threads(num_of_threads/2)
    for(int i=1;i<=slice_total;i++)
    {
        auto covariance_matrix=covariance_matrix_slice(groundstate,i,startpoing);
        entropy=entropy_calculation(covariance_matrix);
        #pragma omp critical
        {
        fp<<i<<"\t"<<entropy<<endl;
        }
    }
    fp.close();
}
void verify_entropy_simple_majoranaGamma(int slice_total, MatrixXcd & Gamma)
{
    double entropy;
    ofstream fp;
    fp.open("./plot_data/" + file_name + "EntropyD" + to_string(D) + "L" + to_string(L) + "alpha" + to_string(alpha) + "uniform" + to_string(onsite_random) + ".txt");
    int * startpoint;
    startpoint=new int[D]{};
    memset(startpoint,0,sizeof(int)*D);
    #pragma omp parallel for private(entropy) schedule(dynamic) num_threads(num_of_threads/2)
    for(int i=1;i<=slice_total;i++)
    {
        auto gamma_matrix=covariance_matrix_slice(Gamma,i*2,startpoint);
        entropy=entropy_calculation_from_MajoranaGamma(gamma_matrix);
        #pragma omp critical
        {
            fp<<i<<"\t"<<entropy<<endl;
        }
    }
    fp.close();
}
void Majorana_hamiltonian(int n)
{
    MatrixXcd cova=MatrixXcd::Zero(Nsites,Nsites);
    L=L/2;
    generating_majoranamatrix(cova,n);
    L=L*2;
    cova=complex<double>(0,1)*cova;
    cova=(cova+cova.adjoint())/2.0;
    file_name="Majorana"+to_string(n);
    MatrixXcd H;
    H=-cova*2.0+MatrixXcd::Identity(Nsites,Nsites);
    verify_entropy_simpleH(Nsites/2,H);
}
void Majorana_gamma(int n)
{
    MatrixXcd gamma=MatrixXcd::Zero(Nsites,Nsites);
    L=L/2;
    generating_majoranamatrix(gamma,n);
    gamma = normalized_gamma_matrix(gamma);
    cout << gamma << endl;
    L=L*2;
    file_name="Majoranagamma"+to_string(n);
    verify_entropy_simple_majoranaGamma(Nsites/4,gamma);
    cout<<"The entropy of the total system is "<<entropy_calculation_from_MajoranaGamma(gamma)<<endl;
}

void verify_entropy(int slice_total)
{
    int slice;
    double entropy;
    ofstream fp;
    int repetition_time = 50;
    fp.open("./plot_data/" + file_name + "EntropyD" + to_string(D) + "L" + to_string(L) + "alpha" + to_string(alpha) + "onsiterandom" + to_string(onsite_random) + ".txt");
    MatrixXd Hamiltonian = MatrixXd::Zero(Nsites, Nsites);
    initialization_node();
    initialization_Hamiltonian_double_2D(Hamiltonian);
    auto groundstate = ground_state(Hamiltonian);
    int starting_point[3] = {0, 0, 0};
    for (int i = 1; i <= slice_total; i++)
    {
        slice = i;
        entropy = 0;
        for (int k = 0; k <= repetition_time - 1; k++)
        {
            for (int j = 0; j < D; j++)
            {
                starting_point[j] = rand() % L;
            }
            auto covariance_matrix = covariance_matrix_slice(groundstate, slice, starting_point);
            entropy += entropy_calculation(covariance_matrix);
        }
        entropy = entropy / repetition_time;
        fp << slice << "\t" << entropy << endl;
    }
    fp.close();
}
void verify_entropy_complex(int slice_total)
{
    int slice;
    double entropy;
    ofstream fp;
    int repetition_time = 50;
    fp.open("./plot_data/" + file_name + "EntropyD" + to_string(D) + "L" + to_string(L) + "alpha" + to_string(alpha) + "onsiterandom" + to_string(onsite_random) + ".txt");
    MatrixXcd Hamiltonian = MatrixXcd::Zero(Nsites, Nsites);
    initialization_node();
    initialization_Hamiltonian_complex(Hamiltonian);
    auto groundstate = ground_state(Hamiltonian);
    int starting_point[3] = {0, 0, 0};
    for (int i = 1; i <= slice_total; i++)
    {
        slice = i;
        entropy = 0;
        for (int k = 0; k <= repetition_time - 1; k++)
        {
            for (int j = 0; j < D; j++)
            {
                starting_point[j] = rand() % L;
            }
            auto covariance_matrix = covariance_matrix_slice(groundstate, slice, starting_point);
            entropy += entropy_calculation(covariance_matrix);
        }
        entropy = entropy / repetition_time;
        fp << slice << "\t" << entropy << endl;
    }
    fp.close();
}
void calculate_number_of_particlepair_different_length(void)
{
    vector<int> system_index;
    int length_per_axis=L;
    int product=1;
    for(int i=0;i<D;i++)
    {
        product*=length_per_axis;
    }
    Nsites=product; 
    initialization_node();
    int lower_index=length_per_axis/2-4;
    int upper_index=length_per_axis/2+5;
    int length_of_system=upper_index-lower_index+1;
    product=1;
    for(int i=0;i<D;i++)
    {
        product*=length_of_system;
    }
    int * temp_array;
    temp_array = new int[D];
    for(int i=0;i<product;i++)
    {
        int i_copy=i;
        for(int j=0;j<D;j++)
        {
            temp_array[j]=i_copy%length_of_system+lower_index;
            i_copy/=length_of_system;
        }
        int temp_index=0;
        product=1;
        for(int j=0;j<D;j++)
        {
            temp_index+=temp_array[j]*product;
            product*=L;
        }
        system_index.push_back(temp_index);
    }
    delete[] temp_array;
    int max_distance=0;
    for(int i=0;i<D;i++)
    {
        max_distance+=length_per_axis;
    }
    int * recording_array=new int[max_distance+1];
    memset(recording_array,0,sizeof(int)*(max_distance+1));
    for(int i=0;i<product;i++)
    {
        for(auto j:system_index) 
        {
            int distance=0;
            for(int k=0;k<D;k++)
            {
                distance+=abs(nodes[i].position[k]-nodes[j].position[k]);
            }
            recording_array[distance]++;
        }
    }
    ofstream fp;
    fp.open("./plot_data/number_of_particlepair_different_length.txt");
    for(int i=0;i<=max_distance;i++)
    {
        if(recording_array[i]!=0)
        fp<<i<<"\t"<<recording_array[i]<<endl;
    }
    fp.close();
}
int main(int argc, char **argv)
{/*
    cout << "the dimension is:\n";
    cin >> D;
    cout << "the L is:\n";
    cin >> L;*/
    D=1;
    L=1000;
    int slice;
   // calculate_number_of_particlepair_different_length();
  /*  cout << "the alpha is:\n";
    cin >> alpha;
    cout << "the modified name is:\n";
    cin >> file_name;
    cout << "the slice is:\n";
    cin >> slice;
    cout << "Onsite random:\n";
    cin >> onsite_random;*/
    alpha=1.1;
    file_name="SSH_modified";
    slice=1000;
    onsite_random=0;
    Nsites = pow(L, D);
    x_floquet = 0.15;
    lambda_floquet = 1/2.0*asinh(cos(2*x_floquet)*cos(2*x_floquet)/sin(2*x_floquet))+0.5; 
    //lambda_floquet = 0.8;
    int floquet_scaling = atoi(argv[1]);
    verify_covariance_scaling(floquet_scaling);
    //verify_entropy_simple(slice);
    num_of_threads = omp_get_max_threads();
    Majorana_gamma(floquet_scaling);
    Majorana_hamiltonian(floquet_scaling);
}