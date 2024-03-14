#include <stdio.h>
#include <math.h>

double f(double x, double xp, double xq, double yq)//定义用于求解反设点的函数，xp,xq,yq分别为xP,xQ,yQ
{
    double y;
    y = (xp * x - 1) / sqrt(xp * xp - 2 * xp * x + 1) - (xq * x + yq * sqrt(1 - x * x) - 1) / sqrt(xq * xq - 2 * xq * x + yq * yq - 2 * yq * sqrt(1 - x * x) + 1);
    return y;
}

double bisect_f(double xp, double xq, double yq, double a, double b, double eps)//定义二分法求解f根的函数, a,b为初态区间，eps为精度限制
{
    double x;
    if (f(a, xp, xq, yq) * f(b, xp, xq, yq) > 0)
    {
        printf("不可用二分法求解\n");
            return 2; //在反射点横坐标可能取值之外，用于辨别二分法失效的情形
    }
    else if (f(a, xp, xq, yq) == 0) //端点即零点
        return a;
    else if (f(b, xp, xq, yq) == 0) //端点即零点
        return b;
    else
    {
        while (fabs(a - b) > eps)
        {
            x = (a + b) / 2;
            if (f(x, xp, xq, yq) * f(a, xp, xq, yq) > 0)
                a = x;
            else if (f(x, xp, xq, yq) * f(b, xp, xq, yq) > 0)
                b = x;
            else
                break;
        }
        return x;
    }
}

inline double tangent_point(double xp)//返回切点横坐标
{
    double xt;
}

int main()
{
    double xp=0, xq=0, yq=0, xt=0, yt=0, xr=0, yr=0, a=-1.0, b=0.0, eps=pow(10,-6);
    printf("观察点横坐标xp=");
    scanf("%lf", &xp);
    while (xp!=100)
    {
        printf("物点横坐标xq=");
        scanf("%lf", &xq);
        printf("物点纵坐标yq=");
        scanf("%lf", &yq);
        getchar();
        xt = bisect_f(xp, xq, yq, a, b, eps);
        if (xt != 2)
        {
            yt = sqrt(1 - xt * xt);
            xr = (xt * xt - xt*(xq * yt + xp * yq - xp * yt) + xq * xp * yt) / (xp * yt);
            yr = yq - (xq * yt) / xt + (xr * yt) / xt;
            printf("反射点T坐标(xt,yt)=(%.6lf,%.6lf)\n", xt, yt);
            printf("像点R的坐标(xr,yr)=(%.6lf,%.6lf)\n", xr, yr);
            printf("\n");
        }
        printf("物点横坐标xp=");
        scanf("%lf", &xp);
    }
    getchar();
    return 0;
}