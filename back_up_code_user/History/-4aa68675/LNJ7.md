# <center> 第一次上机作业实验报告</center>

<center> <b> 戴梦佳 PB20511879 </b> </center>

## 题目

输入观察点坐标$P(x_P,y_P)$和物点坐标$Q(x_Q,y_Q)$。输出反射点$T(x_T,y_T)$，像点$R(x_R,y_R)$。



## 算法

### 求解反射点

<img src="./屏幕截图 2023-03-26 210516.png" alt="屏幕截图 2023-03-26 210516" style="zoom: 33%;" />

使用向量点乘表示夹角余弦，列出方程反射角等于入射角，求解零点从而求解反射点。具体推导如下：
$$
\begin{aligned}
&由几何关系知，TV=OT=(x_T,y_T), TP=(x_P-x_T,-y_T), TQ=(x_Q-x_T,y_Q-y_T)\\
& cos( \angle PTV)=\frac{TP\cdot TV}{|TP||TV|}=\frac{x_P\cdot x_T-1}{\sqrt{(x_P-x_T)^2+y_T^2}}\\
& cos( \angle QTV)=\frac{TQ\cdot TV}{|TQ||TV|}=\frac{x_Q\cdot x_T+y_Q\cdot y_T-1}{\sqrt{(x_P-x_T)^2+(y_Q-y_T)^2}}\\
&令cos( \angle PTV)-cos( \angle QTV)=0并代入x_T^2+y_T^2=1化简\\
& 得方程f(x_T)=\frac{x_P\cdot x_T-1}{\sqrt{(x_P-x_T)^2+(1-x_T^2)}}-\frac{x_Q\cdot x_T+y_Q\cdot \sqrt{1-x_T^2} -1}{\sqrt{(x_P-x_T)^2+(y_Q -\sqrt{1-x_T^2})^2}}\\
\end{aligned}
$$
使用二分法求解$f(x)$的根x。这里有一点很重要，用二分法寻根时，我们默认了 $\angle PTV$和$\angle QTV$都是单调变化的，但实际上，这并不总是成立。容易看出，只有当$\angle POT \in [0, \min{(\angle POQ,\angle POG)} ]$时， 上述单调变化的假设才成立。而更近一步的光学分析表明，如果
设G点的横坐标为$x_G$,连线OQ与圆相交点的横坐标为$x_I$。由基本的几何关系我们可以得出, $x_G=1/x_P$, $x_I=-\frac{1}{\sqrt{1+y_Q^2/X_Q^2}}$。
初态区间$[-1, 0]$, 取$\epsilon=10^{-6}$, 输出结果保留小数点后六位。



### 求解像点
$$
\begin{aligned}
&由几何关系知，R在PT延长线上\\
&\frac{y_R-y_T}{x_R-x_T}=\frac{y_T-y_P}{x_T-x_P}\\
&且QR//OT，即\\
&\frac{y_Q-y_R}{x_Q-x_R}=\frac{y_T}{x_T}\\
&联立解得:\\
&x_R=\frac{x_T^2(y_Q-y_P)+x^T(x_Py_T-x_Py_Q-x_Qy_T)+x_Px_Qy_T}{x_Py_T-x_Ty_P}\\
&y_R=y_Q-\frac{y_T}{x_T}(x_Q-x_R)\\
将T坐标代入即可
\end{aligned}
$$

### 具体实现

编写函数使用二分法求解$f(x)$，将根代入$x_R,y_R$表达式。为了运行方便，使用while连续输入坐标，输入q推出循环。



## 运行结果

