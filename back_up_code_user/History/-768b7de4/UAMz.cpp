#include <stdio.h>
#include <math.h>

double f(double x, double xp, double xq, double yq)//����������ⷴ���ĺ�����xp,xq,yq�ֱ�ΪxP,xQ,yQ
{
    double y;
    y = (xp * x - 1) / sqrt(xp * xp - 2 * xp * x + 1) - (xq * x + yq * sqrt(1 - x * x) - 1) / sqrt(xq * xq - 2 * xq * x + yq * yq - 2 * yq * sqrt(1 - x * x) + 1);
    return y;
}

double bisect_f(double xp, double xq, double yq, double a, double b, double eps)//������ַ����f���ĺ���, a,bΪ��̬���䣬epsΪ��������
{
    double x;
    if (f(a, xp, xq, yq) * f(b, xp, xq, yq) > 0)
    {
        printf("�����ö��ַ����\n");
            return 2; //�ڷ������������ȡֵ֮�⣬���ڱ����ַ�ʧЧ������
    }
    else if (f(a, xp, xq, yq) == 0) //�˵㼴���
        return a;
    else if (f(b, xp, xq, yq) == 0) //�˵㼴���
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

int main()
{
    double xp=0, xq=0, yq=0, xt=0, yt=0, xr=0, yr=0, a=-1.0, b=0.0, eps=pow(10,-6);
    printf("�۲�������xp=");
    scanf_s("%lf", &xp);
    while (xp!=100)
    {
        printf("��������xq=");
        scanf("%lf", &xq);
        printf("���������yq=");
        scanf("%lf", &yq);
        getchar();
        xt = bisect_f(xp, xq, yq, a, b, eps);
        if (xt != 2)
        {
            yt = sqrt(1 - xt * xt);
            xr = (xt * xt - xt*(xq * yt + xp * yq - xp * yt) + xq * xp * yt) / (xp * yt);
            yr = yq - (xq * yt) / xt + (xr * yt) / xt;
            printf("�����T����(xt,yt)=(%.6lf,%.6lf)\n", xt, yt);
            printf("���R������(xr,yr)=(%.6lf,%.6lf)\n", xr, yr);
            printf("\n");
        }
        printf("��������xp=");
        scanf_s("%lf", &xp);
    }
    getchar();
    return 0;
}