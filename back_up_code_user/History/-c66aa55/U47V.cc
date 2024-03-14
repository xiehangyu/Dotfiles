#include <cstdio>
#include <random>

#include <gem5/m5ops.h>

void daxpy(double *X, double *Y, double alpha, const int N)
{
    for (int i = 0; i < N; i++)
    {
        Y[i] = alpha * X[i] + Y[i];
    }
}

void daxsbxpxy(double *X, double *Y, double alpha, double beta, const int N)
{
    for (int i = 0; i < N; i++)
    {
        Y[i] = alpha * X[i] * X[i] + beta * X[i] + X[i] * Y[i];
    }
}

void stencil(double *Y, double alpha, const int N)
{
    for (int i = 1; i < N-1; i++)
    {
        Y[i] = alpha * Y[i-1] + Y[i] + alpha * Y[i+1];
    }
}


void daxpy_unroll(double *X, double *Y, double alpha, const int N)
{
    int i;
    int upperbound=N/4*4;//Here we unroll the loop by 4,because most modern machine can execute 4 float point operations at one instruction
    for(i=0;i<upperbound;i+=4)
    {
        Y[i] = alpha * X[i] + Y[i];
        Y[i+1] = alpha * X[i+1] + Y[i+1];
        Y[i+2] = alpha * X[i+2] + Y[i+2];
        Y[i+3] = alpha * X[i+3] + Y[i+3];
    }
    for(;i<N;i++)//Here we deal with the rest of the loop
    {
        Y[i] = alpha * X[i] + Y[i];
    }
}


void daxsbxpxy_unroll(double *X, double *Y, double alpha, double beta, const int N)
{
    int i;
    int upperbound=N/4*4;//Here we unroll the loop by 4,because most modern machine can execute 4 float point operations at one instruction
    for (i = 0; i < upperbound; i+=4)
    {
        Y[i] = X[i] * (X[i]*alpha + beta  + Y[i]);
        Y[i+1] = X[i+1]* (X[i+1]*alpha + beta  +  Y[i+1]);
        Y[i+2] = X[i+2] *(X[i+2]*alpha + beta  +  Y[i+2]);
        Y[i+3] = X[i+3] *(X[i+3]*alpha + beta  +  Y[i+3]);
    }
    for(;i<N;i++)//Here we deal with the rest of the loop
    {
        Y[i] = X[i] *(X[i]*alpha + beta  +  Y[i]);
    }
}


void stencil_unroll(double *Y, double alpha, const int N)
{
//We also assume that the modern machine can execute 4 float point operations at one instruction

    int upperbound=(N-1)/4*4;
    int i;
    for(i=1;i<upperbound;i+=4)
    {
        Y[i] = alpha * (Y[i-1]+Y[i+1]) + Y[i] ;
        Y[i+1]=alpha * (Y[i]+Y[i+2]) + Y[i+1];
        Y[i+2]=alpha * (Y[i+1]+Y[i+3]) + Y[i+2];
        Y[i+3]=alpha * (Y[i+2]+Y[i+4]) + Y[i+3];
    }
    for(;i<N-1;i+=1)
    {
        Y[i] = alpha * (Y[i-1]+Y[i+1]) + Y[i] ;
    }
}

void stencil_softpipe(double *Y, double alpha, const int N)
{
    double new3,new2,new1,new4,new5,new6,new7,new8;
    new1=Y[5];
    int upperbound=(N-1)/4*4;
    int i;
    for(i=1;i<upperbound-4;i+=4)
    {
        new2=Y[i+8];
        Y[i]=alpha*(Y[i-1]+new2)+Y[i];
        Y[i+1]=alpha*(Y[i]+new3)+new2;
        Y[i+2]=alpha*(Y[i+1]+new4)+new3;
        Y[i+3]=alpha*(Y[i+2]+new1)+new4;
        new1=new8;
        new2=new5;
        new3=new6;
        new4=new7;       
    }
    for(;i<N-1;i+=1)
    {
        Y[i] = alpha * (Y[i-1]+Y[i+1]) + Y[i] ;
    }
}


void simpletestfunction(double *X, double *Y, double alpha,double beta, const int N)
{
    double *Xorigin = new double[N], *Yorigin = new double [N], *Xtest = new double[N], *Ytest = new double[N];
    double totaldifference;
    printf("Test for unroll function daxpy\n");
    for (int i = 0; i < N; i++)
    {
        Xorigin[i] = X[i];
        Yorigin[i] = Y[i];
        Xtest[i] = X[i];
        Ytest[i] = Y[i];
    }
    daxpy(Xorigin,Yorigin,alpha,N);
    daxpy_unroll(Xtest,Ytest,alpha,N);
    totaldifference=0;
    for (int i = 0; i < N; i++)
    {
        totaldifference+=(Yorigin[i]-Ytest[i])* (Yorigin[i]-Ytest[i]);
    }
    printf("The total variance between the roll one and unroll one is %lf\n",totaldifference);
    printf("Test for unroll function daxsbxpxy\n");
    for (int i = 0; i < N; i++)
    {
        Xorigin[i] = X[i];
        Yorigin[i] = Y[i];
        Xtest[i] = X[i];
        Ytest[i] = Y[i];
    }
    daxsbxpxy(Xorigin,Yorigin,alpha,beta,N);
    daxsbxpxy_unroll(Xtest,Ytest,alpha,beta,N);
    totaldifference=0;
    for (int i = 0; i < N; i++)
    {
        totaldifference+=(Yorigin[i]-Ytest[i])* (Yorigin[i]-Ytest[i]);
    }
    printf("The total variance between the roll one and unroll one is %lf\n",totaldifference);
    printf("Test for unroll function stencil\n");
    for (int i = 0; i < N; i++)
    {
        Xorigin[i] = X[i];
        Yorigin[i] = Y[i];
        Xtest[i] = X[i];
        Ytest[i] = Y[i];
    }
    stencil(Yorigin,alpha,N);
    stencil_unroll(Ytest,alpha,N);
    totaldifference=0;
    for (int i = 0; i < N; i++)
    {
        totaldifference+=(Yorigin[i]-Ytest[i])* (Yorigin[i]-Ytest[i]);
    }
    printf("The total variance between the roll one and unroll one is %lf\n",totaldifference);
    printf("Test for softpipe function stencil\n");
    for (int i = 0; i < N; i++)
    {
        Xorigin[i] = X[i];
        Yorigin[i] = Y[i];
        Xtest[i] = X[i];
        Ytest[i] = Y[i];
    }
    stencil(Yorigin,alpha,N);
    stencil_softpipe(Ytest,alpha,N);
    totaldifference=0;
    for (int i = 0; i < N; i++)
    {
        totaldifference+=(Yorigin[i]-Ytest[i])* (Yorigin[i]-Ytest[i]);
    }
    printf("The total variance between the roll one and softpiped one is %lf\n",totaldifference);
    printf("Test for unroll function daxsbxpxy_unroll2\n");
    for (int i = 0; i < N; i++)
    {
        Xorigin[i] = X[i];
        Yorigin[i] = Y[i];
        Xtest[i] = X[i];
        Ytest[i] = Y[i];
    }
    daxsbxpxy(Xorigin,Yorigin,alpha,beta,N);
}
int main()
{
    const int N = 10000;
    double *X = new double[N], *Y = new double[N], alpha = 0.5, beta = 0.1;

    //std::random_device rd;
    std::mt19937 gen(0);
    std::uniform_real_distribution<> dis(1, 2);
    for (int i = 0; i < N; ++i)
    {
        X[i] = dis(gen);
        Y[i] = dis(gen);
    }
    simpletestfunction(X,Y,alpha,beta,N);
    m5_dump_reset_stats(0, 0);
    daxpy(X, Y, alpha, N);
    m5_dump_reset_stats(0, 0);
    daxpy_unroll(X, Y, alpha, N);
    m5_dump_reset_stats(0, 0);
    daxsbxpxy(X, Y, alpha, beta, N);
    m5_dump_reset_stats(0, 0);
    daxsbxpxy_unroll(X, Y, alpha, beta, N);
    m5_dump_reset_stats(0, 0);
    stencil(Y, alpha, N);
    m5_dump_reset_stats(0, 0);
    stencil_unroll(Y, alpha, N);
    m5_dump_reset_stats(0, 0);
    stencil_softpipe(Y, alpha, N);
    m5_dump_reset_stats(0, 0);

    double sum = 0;
    for (int i = 0; i < N; ++i)
    {
        sum += Y[i];
    }
    printf("%lf\n", sum);
    return 0;
}
