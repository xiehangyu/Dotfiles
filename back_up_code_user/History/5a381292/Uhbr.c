//vscode��Ҫʹ��GBK����
#include<stdio.h>
#include<math.h>
double F(double x,double a,double b,double c)//������ַ�����������ĺ���
    {
        double y;
        y=(a*x-1)/sqrt(a*a-2*a*x+1)-(b*x+c*sqrt(1-x*x)-1)/sqrt(b*b-2*b*x+c*c-2*c*sqrt(1-x*x)+1);
        return y;
    }
double G(double x,double a)//������ַ����е�����ĺ���
{
    double y;
    y=1-a*x;
    return y;
}
int main()//������
{   
    double F(double x,double a,double b,double c);//�Զ��庯�����������ַ����������꣩
    double G(double x,double a);//�Զ��庯�����������ַ����е����꣩
    double x_p,x_q,y_q,x_t,y_t,x_r,y_r,a,b,c,d,x,y,epsilon;
    int i;//��������
    a=-1,b=0;//�����������ʼ��������
    epsilon=0.0000000000001;
    printf("������۲��P�ĺ�����x_p=");
    while(scanf("%lf",&x_p)!=EOF)//�������㣬��Ctrl+z+����Enter�˳�����
    {
        printf("���������Q�ĺ�����x_q=");
        scanf("%lf",&x_q);
        printf("���������Q��������y_q=");
        scanf("%lf",&y_q);
        getchar();
        c=G(a,x_p);
        d=G(b,x_p);
        //printf("c=%lf\nd=%lf\n",c,d);
        if(c*d<0)//���ַ����е�����
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
            if(sqrt(1-x*x)<y_q)//ȷ�����ַ�����
                b=x;
            else
                b=-sqrt(1-y_q*y_q);
            a=-1;
            //printf("b=%lf\n",b);
        }
        else
        {
            printf("�޷�ʹ�ö��ַ�Ѱ���е㣬�������ֵ������߸����㷨\n");
        }
        d=F(b,x_p,x_q,y_q);
        c=F(a,x_p,x_q,y_q);
        i=0;
        if(c*d<0)//���ַ����������
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
            x_r=(x_q*(x-x_p)*y-x_p*x*y-y_q*x*(x-x_p))/(-x_p*y);//���������
            y_r=y*(x_r-x_q)/x+y_q;
            printf("������������������=%d\n",i);//������
            printf("�����T������(x_t,y_t)=(%.10lf,%.10lf)\n",x,y);
            printf("���R������(x_r,y_r)=(%.10lf,%.10lf)\n",x_r,y_r);
        }
        else
        {
            printf("�޷�ʹ�ö��ַ����������ֵ������߸����㷨\n");
        }
        printf("\n������۲��P�ĺ�����x_p=");//�ٴμ���
        a=-1,b=0;//�ع���ַ���ʼ����
    }
    getchar();
    return 0;
}