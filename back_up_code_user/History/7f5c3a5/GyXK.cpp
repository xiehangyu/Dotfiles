#include <iostream>
#include <fstream>
#include <physics.hpp>
static int N=200;
using namespace std;
int the_multiplicity_of_count_two_sum(int NA, int sum)
{
    if ((sum<=1)||(sum>=2*NA+1))
    {
        return 0;
    }
    if (sum<=NA+1)
    {
        return sum-1;
    }
    if (sum>NA+1)
    {
        return 2*NA-sum+1;
    }
}
int the_multplicity_of_count_two_minus(int NA,int minus)
{
    if ((minus<=-NA)||(minus>=NA))
    {
        return 0;
    }
    return NA-abs(minus);
}
double Sin_over_sin(int NA,int index)
{
    if (index%(2*N)==0)
    {
        return NA;
    }
    if (index%(2*N)==N)
    {
        if (NA%2==0)
            return -NA;
        else
            return NA;
    }
    double k=index*2*M_PI/(N+0.0);
    return sin(NA*k/2.0)/sin(k/2.0);
}

void sum_over_k(int NA,ofstream & fp)
{
    int i,j,l,m;
    double sum=0;
    for(i=-N/2;i<=N/2-1;i++)
    {
        for(l=-N/2;l<=N/2-1;l++)
      //  for(m=-N/2;m<=N/2-1;m++)
        for(j=-N/2;j<=N/2-1;j++)
        {
            double temp;
            temp=Sin_over_sin(NA,i)*Sin_over_sin(NA,j)*Sin_over_sin(NA,i+j+l)*Sin_over_sin(NA,l);
            //temp=Sin_over_sin(NA,i)*Sin_over_sin(NA,i)*Sin_over_sin(NA,j)*Sin_over_sin(NA,j)*Sin_over_sin(NA,l)*Sin_over_sin(NA,i+j)*Sin_over_sin(NA,i+j+l);
        //    temp=Sin_over_sin(NA,i);
            sum += fabs(temp)*(fabs(i)+fabs(j)+fabs(l))/N;
        }
    }
    fp<<NA<<"\t"<<sum<<endl;
}


void diagram_1(int NA,ofstream & fp)
{
    int s,m;
    double k,sum=0;
    for(s=0;s<=N-1;s++)
    {
        k=2*M_PI*s/(N+0.0);
        for(m=-NA+1;m<=NA-1;m++)
        {
            sum+=cos(k*m)*the_multplicity_of_count_two_minus(NA,m);
        }
    }
    fp<<NA<<"\t"<<sum<<endl;
}

void different_NA(int NA, ofstream & fp)
{
    int count=0;
    int i,j,k,l,m,n;
    for(i=1;i<=NA;i++)
    {
        for(j=1;j<=NA;j++)
        {
            for(k=1;k<=NA;k++)
            {
                for(l=1;l<=NA;l++)
                {
                    for(m=1;m<=NA;m++)
                    {
                        for(n=1;n<=NA;n++)
                            if(((i+j-k-l)%N==0) &&((i+j-m-n)%N==0))
                            {
                                count++;
                            }
                    }
                }
            }
        }
    }
    fp<<NA<<"\t"<<count<<endl;
}



int main()
{
    int max_NA;
    cout<<"input the N you want\n";
    cin>>N;
    cout<<"input max NA\n";
    cin>>max_NA;
    ofstream fp;
    fp.open("./not_important_data/an_idea_trial.txt");
    for(int NA=1;NA<=max_NA;NA++)
    {
        cout<<"Now NA="<<NA<<endl;
        sum_over_k(NA,fp);
    }
}