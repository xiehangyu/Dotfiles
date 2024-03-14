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
using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
static const std::complex<double> I(0,1);
static int D=2;
static int L=15;
static double alpha=2;
struct node{
    int index;
    int * position;
    node(int i)
    {
        index=i;
        position=new int[D];
        for(int j=0;j<D;j++)
        {
            position[j]=i%L;
            i=i/L;
        }
    }
} 
static node * nodes;
static int Nsites=pow(L,D);
static MatrixXd Hamiltonian;
double distance(node node1,node node2)
{
    double distance=0;
    for(int i=0;i<D;i++)
    {
        distance+=pow(node1.position[i]-node2.position[i],2);
    }
    return sqrt(distance);
}
void initialization()
{
    nodes = new node[Nsites];
    for(int i=0;i<Nsites;i++)
    {
        nodes[i]=node(i);
    }
    for(int i=0;i<Nsites;i++)
    {
        for(int j=i+1;j<Nsites;j++)
        {
            Hamiltonian(i,j)=1.0/pow(distance(nodes[i],nodes[j]),alpha);
            Hamiltonian(j,i)=Hamiltonian(i,j);
        }
    }
}