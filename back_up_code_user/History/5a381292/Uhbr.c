//vscode中要使用GBK编码
#include<stdio.h>
#include<math.h>
double F(double x,double a,double b,double c)//定义二分法求反射点横坐标的函数
    {
        double y;
        y=(a*x-1)/sqrt(a*a-2*a*x+1)-(b*x+c*sqrt(1-x*x)-1)/sqrt(b*b-2*b*x+c*c-2*c*sqrt(1-x*x)+1);
        return y;
    }
double G(double x,double a)//定义二分法求切点坐标的函数
{
    double y;
    y=1-a*x;
    return y;
}
int main()//主函数
{   
    double F(double x,double a,double b,double c);//自定义函数声明（二分法求反射点横坐标）
    double G(double x,double a);//自定义函数声明（二分法求切点坐标）
    double x_p,x_q,y_q,x_t,y_t,x_r,y_r,a,b,c,d,x,y,epsilon;
    int i;//迭代次数
    a=-1,b=0;//反射点横坐标初始迭代区间
    epsilon=0.0000000000001;
    printf("请输入观察点P的横坐标x_p=");
    while(scanf("%lf",&x_p)!=EOF)//连续计算，按Ctrl+z+两次Enter退出程序
    {
        printf("请输入物点Q的横坐标x_q=");
        scanf("%lf",&x_q);
        printf("请输入物点Q的纵坐标y_q=");
        scanf("%lf",&y_q);
        getchar();
        c=G(a,x_p);
        d=G(b,x_p);
        //printf("c=%lf\nd=%lf\n",c,d);
        if(c*d<0)//二分法求切点坐标
        {
            while(fabs(a-b)>epsilon)
            {
                x=(a+b)/2;
                //printf("x=%lf\n",x);
                if(G(x,x_p)<0)
                {
                    a=x;
                }
                if(G(x,x_p)>0)
                {
                    b=x;
                }
                else
                {
                    b=x;
                    break;
                }
            }
            if(sqrt(1-x*x)<y_q)//确定二分法区间
                b=x;
            else
                b=-sqrt(1-y_q*y_q);
            a=-1;
            //printf("b=%lf\n",b);
        }
        else
        {
            printf("无法使用二分法寻找切点，请更换初值区间或者更换算法\n");
        }
        d=F(b,x_p,x_q,y_q);
        c=F(a,x_p,x_q,y_q);
        i=0;
        if(c*d<0)//二分法求反射点坐标
        {
            while(fabs(a-b)>epsilon)
            {
                x=(a+b)/2;
                //printf("%lf",x);
                if(F(x,x_p,x_q,y_q)*F(a,x_p,x_q,y_q)<0)
                {
                    b=x;
                }
                if(F(x,x_p,x_q,y_q)*F(b,x_p,x_q,y_q)<0)
                {
                    a=x;
                }
                i=i+1;
            }
            y=sqrt(1-x*x);
            x_r=(x_q*(x-x_p)*y-x_p*x*y-y_q*x*(x-x_p))/(-x_p*y);//求像点坐标
            y_r=y*(x_r-x_q)/x+y_q;
            printf("迭代结束，迭代次数=%d\n",i);//输出结果
            printf("反射点T的坐标(x_t,y_t)=(%.10lf,%.10lf)\n",x,y);
            printf("像点R的坐标(x_r,y_r)=(%.10lf,%.10lf)\n",x_r,y_r);
        }
        else
        {
            printf("无法使用二分法，请更换初值区间或者更换算法\n");
        }
        printf("\n请输入观察点P的横坐标x_p=");//再次计算
        a=-1,b=0;//回归二分法初始区间
    }
    getchar();
    return 0;
}