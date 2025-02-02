---
title: "Your Document Title"
author: "Your Name"
output: pdf_document
header-includes:
   - \usepackage{amsmath}
   - \usepackage{hyperref}
---

# 一维单原子链的振动

运动方程：对第$n$个原子可以写为
$$
m\frac{d^{2}u_{n}}{dt^{2}}=\beta(u_{n+1}+u_{n-1}-2u_{n}) \label{21}
$$
其中$u_{n}$是位移。

根据公式\eqref{eq:21}，...

1.  物理模型：

    -   每一个原胞中只含一个原子（单原子）；

    -   只考虑最近邻相互作用；

    -   可以抽象为质量为$m$的小球被弹性系数为$\beta$的无质量弹簧连接起来的弹性链。

2.  运动方程：对第$n$个原子可以写为
\begin{equation}
m\frac{d^{2}u_{n}}{dt^{2}}=\beta(u_{n+1}+u_{n-1}-2u_{n})
\label{eq:motionalequa}
\end{equation} $$
    m\frac{d^{2}u_{n}}{dt^{2}}=\beta(u_{n+1}+u_{n-1}-2u_{n}) \label{eq:motionalequa}
    \end{equation}
    其中$u_{n}$是位移。

3.  运动方程的解：猜测行波解
    $$u_{n}=Ae^{i(\omega t-naq)}\label{eq:caicejie}$$
    其中$A$与$\omega$待定，$q$为一个参数（这个参数被赋予波矢的意义）。

4.  色散关系：即$\omega$与$q$的关系。将式[\[eq:caicejie\]](#eq:caicejie){reference-type="ref"
    reference="eq:caicejie"}代入式[\[eq:motionalequa\]](#eq:motionalequa){reference-type="ref"
    reference="eq:motionalequa"}解得$\omega$与$q$的关系为
    $$\omega=2\sqrt{\frac{\beta}{m}}|\sin\frac{1}{2}aq|\label{eq:sesan}$$

    -   无论是式[\[eq:caicejie\]](#eq:caicejie){reference-type="ref"
        reference="eq:caicejie"}还是式[\[eq:sesan\]](#eq:sesan){reference-type="ref"
        reference="eq:sesan"}，当$q\to q+\frac{2\pi}{a}$后，两式都不变（注意式[\[eq:sesan\]](#eq:sesan){reference-type="ref"
        reference="eq:sesan"}有一个绝对值），所以，运动方程的解关于$q$具有$\frac{2\pi}{a}$的周期性。[因此，有物理意义的$q$的取值可以被限定在第一布里渊区$-\frac{\pi}{a}\leq q<\frac{\pi}{a}$!]{style="color: red"}

    -   进一步假设，体系中总共有$N$个原子，并采用周期性边界条件，即$u_{n}=u_{n+N}$,代入[\[eq:caicejie\]](#eq:caicejie){reference-type="ref"
        reference="eq:caicejie"}得到$e^{-iNaq}=1$，[
        这要求$q$只能在第一布里渊区中取离散的值$q=\frac{2\pi}{Na}m$,
        $-\frac{N}{2}\leq m<\frac{N}{2}$，$q$可以取的值的总数为$N$。]{style="color: red"}

    -   $q$的取值范围在$[-\frac{\pi}{a},\frac{\pi}{a})$,
        一共有$N$个取值，所以$q$的分布密度为 $$\rho(q)=\frac{Na}{2\pi}$$

    -   $q$被称为准动量，注意和动量的区别！($q$的取值离散，也只能在第一布里渊区中取值，第二点相当于说$q$有平移$\frac{2\pi}{a}$后的不变性。)

    -   相速度定义为$v_{p}=\frac{\omega}{q}$,群速度定义为$v_{g}=\frac{d\omega}{dq}$.

    -   在长波近似区域$q\to0$,
        式[\[eq:sesan\]](#eq:sesan){reference-type="ref"
        reference="eq:sesan"}化为线性色散关系。故此色散关系称为声学支。

# 一维双原子链的振动

1.  物理模型：

    -   一个晶胞中有两个原子。晶格常数为$2a$，每个原子间距为$a$。

    -   这两个原子质量分别为$m$和$M$,前者位于奇数格点上($2n+1$)，后者位于偶数格点上($2n$)。

    -   依旧将相互作用假设为弹簧，且所有原子间弹簧的弹性系数都是$\beta$。

2.  运动方程：对质量为$M$的原子，其在偶数格点上，方程为
    $$M\frac{d^{2}u_{2n}}{dt^{2}}=\beta(u_{2n+1}+u_{2n-1}-2u_{2n})\label{eq:Mequas}$$
    对质量为$m$的原子，其在奇数格点上，方程为
    $$m\frac{d^{2}u_{2n+1}}{dt^{2}}=\beta(u_{2n}+u_{2n+2}-2u_{2n+1})\label{eq:mequas}$$

3.  运动方程的猜测解：由于奇数格点和偶数格点上的原子质量不同。所以，猜测他们也具有不同的振动幅度。但假设他们的振动频率和波矢是一样的。将式[\[eq:caicejie\]](#eq:caicejie){reference-type="ref"
    reference="eq:caicejie"}修改为 $$\begin{cases}
    u_{2n} & =Ae^{i(\omega t-2naq)}\\
    u_{2n+1} & =Be^{i(\omega t-(2n+1)aq)}
    \end{cases}$$ 代入式[\[eq:Mequas\]](#eq:Mequas){reference-type="ref"
    reference="eq:Mequas"}和式[\[eq:mequas\]](#eq:mequas){reference-type="ref"
    reference="eq:mequas"}得 $$\begin{cases}
    (M\omega^{2}-2\beta)A+2\beta\cos(aq)B & =0\\
    2\beta\cos(aq)A+(m\omega^{2}-2\beta)B & =0
    \end{cases}\label{eq:determinant}$$

4.  色散关系：式[\[eq:determinant\]](#eq:determinant){reference-type="ref"
    reference="eq:determinant"}有解的条件是行列式为零，由此可以解得$\omega$与$q$的关系。将行列式写出后，是一个关于$\omega^{2}$的二次方程，有两个解，分别记为$\omega_{\pm}$
    $$\omega_{\pm}=\{\frac{\beta}{mM}[(m+M)\pm(m^{2}+M^{2}+2mM\cos2qa)^{\frac{1}{2}}]\}^{\frac{1}{2}}\label{eq:sesan2atom}$$

    -   $q$的取值依然被限定在第一布里渊区内。注意此时，晶格常数为$2a$(一个晶胞中有两个原子），所以，第一布里渊区的取值范围为
        $$-\frac{\pi}{2a}\leq q<\frac{\pi}{2a}$$

        备注：式[\[eq:sesan2atom\]](#eq:sesan2atom){reference-type="ref"
        reference="eq:sesan2atom"}具有平移$q\to q+\frac{\pi}{a}$后的不变形。但式[\[eq:determinant\]](#eq:determinant){reference-type="ref"
        reference="eq:determinant"}似乎不具有。实际上，当$q\to q+\frac{\pi}{a}$后，$u_{2n}\to u_{2n},u_{2n+1}\to-u_{2n+1}$,但这个负号可以被吸收进$B$里面去（至于待定系数究竟是$B$还是$-B$都无所谓。可以理解成这个待定系数的大小最终还是由初始条件去决定的）。所以，这个负号不影响。可以认为，依然具有平移$q\to q+\frac{\pi}{a}$的不变性。

    -   色散关系曲线见下图

        ![image](\string"Figure/Screenshot from 2023-04-13 22-37-49\string".png){width="45%"}

    -   $\omega_{-}$:元胞的质心运动，声学支；$\omega_{+}$:
        元胞中原子的相对运动，光学支。

    -   周期性边界条件：假设体系中一共有$N$个晶胞（$2N$个原子），则$q$同样只能在第一布里渊区中离散取值，且$e^{-2iNaq}=1$,
        即
        $$q=\frac{\pi}{Na}m,\text{-\ensuremath{\frac{N}{2}\leq m<\frac{N}{2}}}$$
        $q$的密度为$\rho=\frac{Na}{\pi}$。

    -   该体系一共有两支格波($\omega_{\pm}$),
        每支格波中$q$有$N$中取值（对应着晶胞总数），两支格波加一起共有$2N$种振动模式(对应着原子总数)。

# 3维情况

若该体系共有$N$个晶胞，每个晶胞中有$n$个原子，那么该体系

-   共有$3n$支格波（或称$3n$支色散关系），$3$支为声学支，$3n-3$支为光学支。

-   每一支的$q$有$N$种取值，$\rho(q)=\frac{V}{(2\pi)^{3}}$。

-   一共有$3nN$中$3nN$种振动模式。
