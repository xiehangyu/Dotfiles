#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
using namespace std;
void DooliteDecomposition(vector<vector<double>> &A)
//Do the Doolite Decomposition for matrix A, here after the decomposition, the L, U matrix is stored in A
{
    //This function is to calculate the determinant
    double temp_sum;
    int n=A.size();//Get the size of the matrix
    for(int k=0;k<=n-1;k++)
    {
        for(int j=k;j<=n-1;j++)
        //Calculate U, and U is stored in the upper triangle of A
        {
            temp_sum=0;
            for(int s=0;s<=k-1;s++)
            {
                temp_sum+=A[k][s]*A[s][j];
            }
            A[k][j]=A[k][j]-temp_sum;
        }
        for(int i=k+1;i<=n-1;i++)
        //Calculate L, and L is stored in the lower triangle of A
        {
            temp_sum=0;
            for(int s=0;s<=k-1;s++)
            {
                temp_sum+=A[i][s]*A[s][k];
            }
            A[i][k]=(A[i][k]-temp_sum)/A[k][k];
        }
    }
}
void MatrixProduct(vector<vector<double>> & A,vector<vector<double>> & B,vector<vector<double>> & C)
{

//Calculate the matrix product of A*B, the result is stored in C
    int n=A.size();
    int m=B[0].size();
    int l=B.size();
    for(int i=0;i<=n-1;i++)
    {
        for(int j=0;j<=m-1;j++)
        {
            C[i][j]=0;
            for(int k=0;k<=l-1;k++)
            {
                C[i][j]+=A[i][k]*B[k][j];
            }
        }
    }

}
double offdiag(vector<vector<double>> & M)
{
    //Calculate the square value of the off-diagonal elements of the matrix M
    int n=M.size();
    double sum=0;
    for(int i=0;i<=n-1;i++)
    {
        for(int j=0;j<=n-1;j++)
        {
            if(i!=j)
            {
                sum+=M[i][j]*M[i][j];
            }
        }
    }
    return sum;
}
void JacobyEig(vector<vector<double>> & M,vector<vector<double>> &  Q,vector<double> & eigenvalue, int filename=0,bool whethertest=false)
{
// Here we use the Jacoby method to calculate the eigenvalues and eigenvectors of a symmetric matrix, and store them in eigenvalue and Q respectively. M is the input
// matrix, Q is the matrix of eigenvectors, and eigenvalue is the vector of eigenvalues.
    double epsilon=1E-6; // The precision of the calculation
    int n=M.size(); // The dimension of the matrix
    Q.resize(n);
    string filenamestr="./data/off_diagonal_square"+to_string(filename)+".txt";
    ofstream gp;
    gp.open(filenamestr);
    for(int i=0;i<=n-1;i++)
    {
        Q[i].resize(n);
        Q[i][i]=1;
    }
    eigenvalue.resize(n);
    auto M_re=M;
    double max=0; //Store the maximum value of the off-diagonal elements
    int row_max=0; //Store the row number of the maximum value
    int col_max=0; //Store the column number of the maximum value
    int step=0;
    double off_diagonal_sum=offdiag(M_re);
    gp<<step<<"\t"<<off_diagonal_sum<<endl;
    if(whethertest==true)
    {
        cout<<"The square of off-diagonal elements during calculations are:"<<endl;
        cout<<"step\t<<value"<<endl;
        cout<<step<<"\t"<<off_diagonal_sum<<endl;
    }
    while(off_diagonal_sum>epsilon)
    {
        max=0;
        for(int i=0;i<=n-1;i++)
        {
            for(int j=0;j<=n-1;j++)
            {
                if(i!=j)
                {
                    if(abs(M_re[i][j])>max)
                    {
                        max=abs(M_re[i][j]);
                        row_max=i;
                        col_max=j;
                    }
                }
            }
        }
        double s=(M_re[col_max][col_max]-M_re[row_max][row_max])/(2*M_re[row_max][col_max]);
        double t=0;
        if(s>=0)
        {
            t=1/(s+sqrt(1+s*s));
        }
        else
        {
            t=1/(s-sqrt(1+s*s));
        }
        double c=1/sqrt(1+t*t);
        double d=c*t;
        double temp1,temp2;//store a_pi,a_qi
        for(int i=0;i<=n-1;i++)
        {
            if(i!=row_max&&i!=col_max)
            {
                temp1=M_re[row_max][i];
                temp2=M_re[col_max][i];
                M_re[row_max][i]=c*temp1-d*temp2;
                M_re[col_max][i]=d*temp1+c*temp2;
                M_re[i][row_max]=M_re[row_max][i];
                M_re[i][col_max]=M_re[col_max][i];
            }
        }
        M_re[row_max][row_max]-=t*M_re[row_max][col_max];
        M_re[col_max][col_max]+=t*M_re[row_max][col_max];
        M_re[row_max][col_max]=0;
        M_re[col_max][row_max]=0;
        for(int i=0;i<=n-1;i++)
        {
            temp1=Q[i][row_max];
            temp2=Q[i][col_max];
            Q[i][row_max]=c*temp1-d*temp2;
            Q[i][col_max]=d*temp1+c*temp2;
        }//Update Q
        off_diagonal_sum=offdiag(M_re);
        gp<<++step<<"\t"<<off_diagonal_sum<<endl;
        if(whethertest==true)
        {
            cout<<step<<"\t"<<off_diagonal_sum<<endl;
        }
    }
    eigenvalue.resize(n);
    for(int i=0;i<=n-1;i++)
    {
        eigenvalue[i]=M_re[i][i];
    }
}

void createMatrix(vector<vector<double>> & M, vector<int> & label)
{
    //read the data of iris from the file
    string filename="./data/iris.txt";
    ifstream fp;
    fp.open(filename);
    int end=150;
    char c;
    M.resize(4);
    label.resize(150);
    for(int i=0;i<=3;i++)
    {
        M[i].resize(150);
    }
    for(int i=0;i<=end-1;i++)
    {
        fp>>M[0][i]>>c>>M[1][i]>>c>>M[2][i]>>c>>M[3][i]>>c>>label[i];
    }
    //Normalize the data
    for(int i=0;i<=3;i++)
    {
        double sum=0;
        for(int j=0;j<=end-1;j++)
        {
            sum+=M[i][j];
        }
        sum=sum/(end+0.0);
        for(int j=0;j<=end-1;j++)
        {
            M[i][j]-=sum;
        }
    }
    fp.close();
}
vector<vector<double>> transpose_matrix(vector<vector<double>> & M)
//Calculate the transpose of the matrix M
{
    int m=M.size();
    int n=M[0].size();
    vector<vector<double>> MT;
    MT.resize(n);
    for(int i=0;i<=n-1;i++)
    {
        MT[i].resize(m);
    }
    for(int i=0;i<=m-1;i++)
    {
        for(int j=0;j<=n-1;j++)
        {
            MT[j][i]=M[i][j];
        }
    }
    return MT;
}
void covariance_matrix_fromM(vector<vector<double>> & M,vector<vector<double>> & cov)
{
    int end=150;
    //Calculate the covariance matrix of the data matrix M
    cov.resize(end);
    for(int i=0;i<=end-1;i++)
    {
        cov[i].resize(end);
    }
    auto MT=transpose_matrix(M);
    MatrixProduct(M,MT,cov);
    for(int i=0;i<=end-1;i++)
    {
        for(int j=0;j<=end-1;j++)
        {
            cov[i][j]=cov[i][j]/(end+0.0);
        }
    }
}
void sorted_eigensystems(vector<double> & eigenvalues, vector<vector<double>> & eigenvectors)
{
    auto temp_matrix=transpose_matrix(eigenvectors);
    int n=eigenvalues.size();
    for(int i=0;i<=n-1;i++)
    {//bubble sort
        for(int j=i+1;j<=n-1;j++)
        {
            if(eigenvalues[i]<eigenvalues[j])
            {
                swap(eigenvalues[i],eigenvalues[j]);
                swap(temp_matrix[i],temp_matrix[j]);
            }
        }
    }
    eigenvectors=transpose_matrix(temp_matrix);
}
double inner_product(vector<double> & v1,vector<double> & v2)
//Calculate the inner product of two vectors
{
    double sum=0;
    int n=v1.size();
    for(int i=0;i<=n-1;i++)
    {
        sum+=v1[i]*v2[i];
    }
    return sum;
}
void PCA()
{
    vector<vector<double>> M;
    vector<int> label;
    createMatrix(M,label);
    vector<vector<double>> cov;
    covariance_matrix_fromM(M,cov);
    cout<<endl;
    cout<<"The covariance matrix is:"<<endl;
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=3;j++)
        {
            cout<<cov[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
    vector<double> eigenvalues;
    vector<vector<double>> eigenvectors;
    JacobyEig(cov,eigenvectors,eigenvalues);
    sorted_eigensystems(eigenvalues,eigenvectors);
    cout<<"The eigenvalues for covariance matrix are:"<<endl;
    for(int i=0;i<=3;i++)
    {
        cout<<eigenvalues[i]<<"\t";
    }
    cout<<endl;
    cout<<"The eigenvectors for covariance matrix are:"<<endl;
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=3;j++)
        {
            cout<<eigenvectors[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
    vector<vector<double>> projection_matrix;
    projection_matrix.resize(2);
    for(int i=0;i<=1;i++)
    {
        projection_matrix[i].resize(4);
    }
    for(int i=0;i<=1;i++)
    {
        for(int j=0;j<=3;j++)
        {
            projection_matrix[i][j]=eigenvectors[j][i];
        }
    }
    auto MT=transpose_matrix(M);
    ofstream fp;
    fp.open("./data/iris_PCA.txt");
    for(int i=0;i<=149;i++)
    {  
        for(int j=0;j<=1;j++)
        {
           // cout<<inner_product(projection_matrix[j],MT[i])<<"\t";
            fp<<inner_product(projection_matrix[j],MT[i])<<"\t";
        }
      //  cout<<label[i]<<endl;
        fp<<label[i]<<endl;
    }
    fp.close();
}

void UVDdecomposition(vector<vector<double>> & M, vector<vector<double>> & U, vector<vector<double>> & D, vector<vector<double>> & V)
{
    auto MT=transpose_matrix(M);
    vector<vector<double>> AAT;
    AAT.resize(M.size());
    for(int i=0;i<=M.size()-1;i++)
    {
        AAT[i].resize(MT[0].size());
    }
    MatrixProduct(M,MT,AAT);
    vector<vector<double>> ATA;
    ATA.resize(MT.size());
    for(int i=0;i<=MT.size()-1;i++)
    {
        ATA[i].resize(M[0].size());
    }
    MatrixProduct(MT,M,ATA);
    vector<double> eigenvalues;
    JacobyEig(AAT,U,eigenvalues,1,true);
    sorted_eigensystems(eigenvalues,U);
    JacobyEig(ATA,V,eigenvalues,2);
    sorted_eigensystems(eigenvalues,V);
    D.resize(M.size());
    for(int i=0;i<=M.size()-1;i++)
    {
        D[i].resize(MT.size());
        for(int j=0;j<=MT.size()-1;j++)
        {
            D[i][j]=0;
        }
    }
    for(int i=0;i<min(M.size(),MT.size());i++)
    {
        D[i][i]=sqrt(eigenvalues[i]);
    }
}
void random_matrix_SVD()
{
    vector<vector<double>> M;
    M.resize(4);
    for(int i=0;i<=3;i++)
    {
        M[i].resize(3);
    }
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=2;j++)
        {
            M[i][j]=(double)rand()/(RAND_MAX+0.0);
        }
    }
    cout<<"The random matrix is:"<<endl;
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=2;j++)
        {
            cout<<M[i][j]<<"\t";
        }
        cout<<endl;
    }
    vector<vector<double>> U;
    vector<vector<double>> D;
    vector<vector<double>> V;
    cout<<endl;
    UVDdecomposition(M,U,D,V);
    cout<<"The SVD decomposition of the random matrix is:"<<endl;
    cout<<"The U matrix is:"<<endl;
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=3;j++)
        {
            cout<<U[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
    cout<<"The Sigma matrix is:"<<endl;
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=2;j++)
        {
            cout<<D[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
    auto UT=transpose_matrix(U);
    auto VT=M;
    MatrixProduct(UT,M,VT);

    cout<<"The transpose of V is:"<<endl;
    for(int i=0;i<=2;i++)
    {
        for(int j=0;j<=2;j++)
        {   
            if(D[i][i]!=0)
                cout<<VT[i][j]/D[i][i]<<"\t";
            else
                cout<<VT[i][j]<<"\t";
        }
        cout<<endl;
    }
}
int main()
{
    cout.precision(3);
    random_matrix_SVD();
    PCA();
}
