#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//#define double float
using namespace std;
void DooliteDecomposition(vector<vector<double>> &A)
//Do the Doolite Decomposition for matrix A, here after the decomposition, the L, U matrix is stored in A
{
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
    //print the decomposition result
    /*
    cout<<"The L matrix is:"<<endl;
    for(int k=0;k<=n-1;k++)
    {
        for(int j=0;j<=k-1;j++)
        {
            cout<<A[k][j]<<"  ";
        }
        cout<<1<<"\n";
    }
    cout<<"The U matrix is:"<<endl;
    for(int k=0;k<=n-1;k++)
    {
        for(int j=0;j<=n-1;j++)
        {
            if(j<k)
            {
                cout<<0<<"  ";
            }
            else
            {
                cout<<A[k][j]<<"  ";
            }
        }
        cout<<"\n";
    }*/
}

void calculate_solution(vector<vector<double>> &A, vector<double> &b, vector<double> &x)
// Solve the equation Ax=b, where A has already been Doolite decomposed
{
    int n=b.size();
    vector<double> y(n);
    //calculate Ly=b
    double temp_sum=0;
    for(int i=0;i<=n-1;i++)
    {
        temp_sum=0;
        for(int j=0;j<=i-1;j++)
        {
            temp_sum+=A[i][j]*y[j];
        }
        y[i]=b[i]-temp_sum;
    }
    cout<<endl;
    //calculate Ux=y
    for(int i=n-1;i>=0;i--)   
    {
        temp_sum=0;
        for(int j=n-1;j>=i+1;j--)
        {
            temp_sum+=A[i][j]*x[j];
        }
        x[i]= (y[i]-temp_sum)/A[i][i];
    }
}

vector<vector<double>> generate_matrixA()
// generate the first matrix A
{
    vector<vector<double>> A(5);
    for(int i=0;i<=4;i++)
    {
        A[i].resize(5);
    }
    for(int i=0;i<=4;i++)
    {
        for(int j=0;j<=4;j++)
        {
            A[i][j]=1.0/(9.0-i-j);
        }
    }
    /*
    cout<<"The first matrix is"<<endl;
    for(int i=0;i<=4;i++)
    {
        for(int j=0;j<=4;j++)
        {
            cout<<A[i][j]<<"  ";
        }
        cout<<"\n";
    }
    */
    return A;
}

vector<vector<double>> generate_matrixB()
//generate the second matrix
{
    vector<vector<double>> B(4);
    B[0]=vector<double>{4.0,-1.0,1.0,3.0};
    B[1]=vector<double>{16.0,-2.0,-2.0,5.0};
    B[2]=vector<double>{16.0,-3.0,-1.0,7.0};
    B[3]=vector<double>{6,-4,2,9};
    /*
    cout<<"The second matrix is"<<endl;
    for(int i=0;i<=3;i++)
    {
        for(int j=0;j<=3;j++)
        {
            cout<<B[i][j]<<"  ";
        }
        cout<<"\n";
    }*/
    return B;
}

void calculate_min_eigenvalue(string filename, vector<vector<double>> &A)
{
    //Calculate the minimum eigenvalue of matrix A and store the middle result in a csv file
    double epislon=1E-5;//set the precision
    ofstream fp;
    fp.open(filename);
    int n=A.size();
    vector<double> y(n);
    vector<double> x(n);
    fp<<"k,u,";
    for(int i=0;i<=n-1;i++)
    {
        fp<<"x"<<i+1<<",";
    }
    for(int i=0;i<=n-2;i++)
    {
        fp<<"y"<<i+1<<",";
    }
    fp<<"y"<<n<<"\n";
    for(int i=0;i<=n-1;i++)
    {
        x[i]=1;
        y[i]=1;
    }
    fp<<0<<","<<1<<",";
    for(int i=0;i<=n-1;i++)
    {
        fp<<x[i]<<",";
    }
    for(int i=0;i<=n-2;i++)
    {
        fp<<y[i]<<",";
    }
    fp<<y[n-1]<<endl;
    double lambda_now=1;
    double lambda_last=0;//lambda is the infinite norm
    int step=0;
    DooliteDecomposition(A);
    do{
        lambda_last=lambda_now;
        calculate_solution(A,y,x);
        step++;
        lambda_now=0;
        for(int i=0;i<=n-1;i++)
        {
            if(fabs(x[i])>fabs(lambda_now))
            {
                //Here we want to keep the minus or positive sign of the eigenvalues
                lambda_now=(x[i]);
            }
        }
        for(int i=0;i<=n-1;i++)
        {
            y[i]=x[i]/lambda_now;
        }
        fp<<step<<","<<lambda_now<<",";
        for(int i=0;i<=n-1;i++)
        {
            fp<<x[i]<<",";
        }
        for(int i=0;i<=n-2;i++)
        {
            fp<<y[i]<<",";
        }
        fp<<y[n-1]<<endl;
    }while(fabs(lambda_now-lambda_last)>epislon);
    //output the eigenvalue and eigenvectors:
    cout<<"The minimum eigenvalue is "<<1.0/lambda_now<<endl;
    cout<<"The corresponding eigenvector is"<<endl;
    for(int i=0;i<=n-1;i++)
    {
        cout<<y[i]<<"\t";
    }
}
void DooliteDecompositionPartialPivot(vector<vector<double>> &A, vector<int> &P)
//Partial Pivot version of Doolite Decomposition, P records the permutation made, L,U is still stored in A
{
    int n=A.size();
    double tempsum=0;
    int max_label;//record the index of pivot
    int max_value;//record the value of pivot
    for(int i=0;i<=n-1;i++)
    {
        //choose the partial pivot
        max_value=0;
        for(int s=i;s<=n-1;s++)
        {
            for(int j=0;j<=i-1;j++)
            {
                tempsum+=A[s][j]*A[j][i];
            }
            if(fabs(A[s][i]-tempsum)>max_value)
            {
                max_value=fabs(A[s][i]-tempsum);
                max_label=s;
            }
        }
        swap(A[i],A[max_label]);//SWAP the rows to dothe partial pivot
        P[i]=max_label;
        //The rest is the same procedure as Doolite Decomposition
        for(int j=i;j<=n-1;j++)
        {
            tempsum=0;
            for(int k=0;k<=i-1;k++)
            {
                tempsum+=A[i][k]*A[k][j];
            }
            A[i][j]=A[i][j]-tempsum;
        }
        for(int j=i+1;j<=n-1;j++)
        {
            tempsum=0;
            for(int k=0;k<=i-1;k++)
            {
                tempsum+=A[j][k]*A[k][i];
            }
            A[j][i]=(A[j][i]-tempsum)/A[i][i];
        }
    }
    /*
    cout<<"The Doolite Decomposition of the matrix is"<<endl;
    for(int i=0;i<=n-1;i++)
    {
        for(int j=0;j<=n-1;j++)
        {
            cout<<A[i][j]<<"  ";
        }
        cout<<"\n";
    }
    cout<<"The permutation matrix is"<<endl;
    for(int i=0;i<=n-1;i++)
    {
        cout<<P[i]<<"  ";
    }
    cout<<"\n";
    */
}
void swap(double &a, double &b)
{
    double temp=a;
    a=b;
    b=temp;
}
void CalculationSolutionDoolitePartialPivot(vector<vector<double>> &A, vector<int> &P, vector<double> &b, vector<double> &x)
//Calculate Ax=b, here the partial pivot Doolite Decomposition is used
{
    //Calculate QLUx=b, equivalent to LUX=Qb
   int n=P.size();
   for(int i=0;i<=n-1;i++)
   {
       swap(b[i],b[P[i]]);
   }
   calculate_solution(A,b,x);
}

void test()//Just a test function to validate the correcteness of the partial pivot Doolite Decomposition
{
    auto A=generate_matrixA();
    auto B=generate_matrixB();
    vector<int> PA(A.size());
    vector<int> PB(B.size());
    vector<double> bA{1,1,1,1,1};
    vector<double> bB{1,1,1,1};
    vector<double> xA(A.size());
    vector<double> xB(B.size());
    DooliteDecompositionPartialPivot(A,PA);
    DooliteDecompositionPartialPivot(B,PB);
    CalculationSolutionDoolitePartialPivot(A,PA,bA,xA);
    CalculationSolutionDoolitePartialPivot(B,PB,bB,xB);
    cout<<"The solution of Ax=b is"<<endl;
    for(int i=0;i<=xA.size()-1;i++)
    {
        cout<<xA[i]<<"  ";
    }
    cout<<"\n";
    cout<<"The solution of Bx=b is"<<endl;
    for(int i=0;i<=xB.size()-1;i++)
    {
        cout<<xB[i]<<"  ";
    }
    cout<<"\n";
}

void calculate_min_eigenvalue_partialpivot(string filename, vector<vector<double>> &A)
//Using partialpivot Doolite Decomposition to calculate the minimum eigenvalue of a matrix
{
    double epislon=1E-5;
    ofstream fp;
    fp.open(filename);
    int n=A.size();
    vector<int> P(n);
    vector<double> y(n);
    vector<double> x(n);
    fp<<"k,u,";
    for(int i=0;i<=n-1;i++)
    {
        fp<<"x"<<i+1<<",";
    }
    for(int i=0;i<=n-2;i++)
    {
        fp<<"y"<<i+1<<",";
    }
    fp<<"y"<<n<<"\n";
    for(int i=0;i<=n-1;i++)
    {
        y[i]=1;
        x[i]=1;
    }
    fp<<0<<","<<1<<",";
    for(int i=0;i<=n-1;i++)
    {
        fp<<x[i]<<",";
    }
    for(int i=0;i<=n-2;i++)
    {
        fp<<y[i]<<",";
    }
    fp<<y[n-1]<<endl;
    double lambda_now=1;
    double lambda_last=0;//lambda is the infinite norm
    int step=0;
    DooliteDecompositionPartialPivot(A,P);
    do{
        lambda_last=lambda_now;
        CalculationSolutionDoolitePartialPivot(A,P,y,x);
        step++;
        lambda_now=0;
        for(int i=0;i<=n-1;i++)
        {
            if(fabs(x[i])>fabs(lambda_now))
            {
                //Here we want to keep the minus or positive signs of the eigenvalues.
                lambda_now=(x[i]);
            }
        }
        for(int i=0;i<=n-1;i++)
        {
            y[i]=x[i]/lambda_now;
        }
        fp<<step<<","<<lambda_now<<",";
        for(int i=0;i<=n-1;i++)
        {
            fp<<x[i]<<",";
        }
        for(int i=0;i<=n-2;i++)
        {
            fp<<y[i]<<",";
        }
        fp<<y[n-1]<<endl;
    }while(fabs(lambda_now-lambda_last)>epislon);
    //output the eigenvalue and eigenvectors:
    cout<<"The minimum eigenvalue is "<<1.0/lambda_now<<endl;
    cout<<"The corresponding eigenvector is"<<endl;
    for(int i=0;i<=n-1;i++)
    {
        cout<<y[i]<<"\t";
    }
}

int main()
{
    auto A1=generate_matrixA();
    auto A2=generate_matrixB();
    calculate_min_eigenvalue("A1Doolite.csv",A1);
    calculate_min_eigenvalue("A2Doolite.csv",A2);
    A1=generate_matrixA();
    A2=generate_matrixB();
    calculate_min_eigenvalue_partialpivot("A1DoolitePartialPivot.csv",A1);
    calculate_min_eigenvalue_partialpivot("A2DoolitePartialPivot.csv",A2);
}
