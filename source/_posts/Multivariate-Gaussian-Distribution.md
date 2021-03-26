---
title: Multivariate Gaussian Distribution
date: 2021-03-24 15:57:31
tags:
mathjax: true
---

# 多元高斯分布

## 引子

先来看一张图：

![图1](./Multivariate-Gaussian-Distribution/1.png)

如果我让你从A、B中选出一个离群点，你觉得谁更合适？从欧式距离来比较，A、B两点到均值（图中的原点）的距离相同，但从直觉上看，B显然更融入群体。这是因为两者在比较时尺度或者说量纲（scalar）没有统一，就像你拿百万面额的津巴布韦币同人民币做比较一样，图中对应的尺度就是维度的方差。当你对每个维度的尺度进行归一化处理后（即让维度上的每个数据除以维度对应的标准差（方差的算数平方根）为什么是平方根？因为我们最终希望修正后的尺度可以简单采用欧式距离，而欧式距离的定义中显然有一个平方根）图就变成下面这样：

![图2](./Multivariate-Gaussian-Distribution/3.jpg)

现在就可以简单地用欧式距离判断谁是离群点了

## 一元高斯分布

如果你学过概率论，那你应该对上面的引子感到熟悉，因为在一元高斯分布中，我们经常会对随机变量$X$进行标准化——$Z=\frac{X-\mu}{\sigma}$，这同之前我们所做的工作相比不能说完全一致，只能说一模一样（XD）。引子就先告一段落，我们很快会再来回顾它。先来复习一下一元高斯分布的有关知识：
$$
\begin{aligned}
p(x) &=\frac{1}{\sigma \sqrt{2 \pi}} \cdot e^{-\frac{1}{2} \cdot\left(\frac{x-\mu}{\sigma}\right)^{2}} \\
1 &=\int_{-\infty}^{+\infty} p(x)dx
\end{aligned}
\tag{1}
$$
令$Z=\frac{X-\mu}{\sigma}$
$$
\begin{aligned}
\because x(z) &=z \cdot \sigma+\mu \\
\therefore p(x(z)) &=\frac{1}{\sigma \sqrt{2 \pi}} \cdot e^{-\frac{1}{2} \cdot(z)^{2}} \\
1 &=\int_{-\infty}^{+\infty} p(x(z)) d x \\
&=\int_{-\infty}^{+\infty} \frac{1}{\sigma \sqrt{2 \pi}} \cdot e^{-\frac{1}{2} \cdot(z)^{2}} d x \\
&=\int_{-\infty}^{+\infty} \frac{1}{\sqrt{2 \pi}} \cdot e^{-\frac{1}{2} \cdot(z)^{2}} dz 
\end{aligned}
\tag{2}
$$
此时我们说随机变量 $Z \sim \mathcal{N}(0,1)$ 服从一元标准高斯分布, 其均值 $\mu=0,$ 方差 $\sigma^{2}=1,$ 其概率密度函数为
$$
p(z)=\frac{1}{\sqrt{2 \pi}} \cdot e^{-\frac{1}{2} \cdot(z)^{2}}
\tag{3}
$$

## 多元高斯分布

### 1.概率密度函数（PDF）

接下来我们将正式探讨多元情况下的高斯分布，从一元到多元，我们的随机变量也进化成了随机向量，我们用$X=(x_1,x_2,...,x_n)$表示，$x_i$是向量的一个分量，是单个随机变量。我们先从最简单的情况，各个随机变量之间相互独立，开始讨论。

假设我们有随机向量 $\vec{Z}=\left[Z_{1}, \cdots, Z_{n}\right]^{\top},$ 其中 $Z_{i} \sim \mathcal{N}(0,1)(i=1, \cdots, n)$ 且 $Z_{i}, Z_{j}(i, j=1, \cdots, n \wedge i \neq j)$ 彼此独立, 即随机向量中的每个随机变量 $Z_{i}$ 都服从标准高斯分布且两两彼此独立. 则由(3)与独立随机变量概率密度函数之间的关系, 我们可得随机向量 $\vec{Z}=\left[Z_{1}, \cdots, Z_{n}\right]^{\top}$ 的联合概率密度函数为
$$
\begin{aligned}
p\left(z_{1}, \cdots, z_{n}\right) &=\prod_{i=1}^{n} \frac{1}{\sqrt{2 \pi}} \cdot e^{-\frac{1}{2} \cdot\left(z_{i}\right)^{2}} \\
&=\frac{1}{(2 \pi)^{\frac{n}{2}}} \cdot e^{-\frac{1}{2} \cdot\left(Z^{\top} Z\right)} \\
1 &=\int_{-\infty}^{+\infty} \cdots \int_{-\infty}^{+\infty} p\left(z_{1}, \cdots, z_{n}\right) d z_{1} \cdots d z_{n}
\end{aligned}
\tag{4}
$$
我们称随机向量 $\vec{Z} \sim \mathcal{N}(\overrightarrow{0}, \mathbf{I})$, 即随机向量服从均值为零向量, 协方差矩阵为单位矩阵的高斯分布. 在这里, 随机向量 $\vec{Z}$ 的协方差矩阵是 $\operatorname{Cov}\left(Z_{i}, Z_{j}\right), i, j=1, \cdots, n$ 组成的矩阵, 即
$$
\begin{aligned}
\left[\operatorname{Cov}\left(Z_{i}, Z_{j}\right)\right]_{n \times n} &=\mathbf{E}\left[(Z-\vec{\mu})(Z-\vec{\mu})^{\top}\right] \\
&=\mathbf{I}
\end{aligned}
\tag{5}
$$
由于随机向量 $\vec{Z} \sim \mathcal{N}(\overrightarrow{0}, \mathbf{I}),$ 所以其协方差矩阵的对角线元素为1, 其余元素为0. 如果我们取常数 $c=p\left(z_{1}, \cdots, z_{n}\right),$ 则可得函数 $p\left(z_{1}, \cdots, z_{n}\right)$ 的等高线为 $c^{\prime}=Z^{\top} Z,$ 当随机向量 $\vec{Z}$ 为二维向量
时, 我们有
$$
c^{\prime}=Z^{\top} \cdot Z=\left(z_{1}-0\right)^{2}+\left(z_{2}-0\right)^{2}
\tag{6}
$$
显然，其等高线是以(0，0)为圆心的同心圆。

![图3](./Multivariate-Gaussian-Distribution/4.jpg)

接下来讨论各随机变量不相互独立的一般情况：

既然我们能够轻松地处理独立的情况，那只要我们能够把一般的情况转化成独立的情况，那么问题是不是就迎刃而解了呢？答案是肯定的，幸运的是，我们有如下定理：

**定理1: 若存在随机向量 $\vec{X} \sim \mathcal{N}(\vec{\mu}, \Sigma),$ 其中 $\vec{\mu} \in R^{n}$ 为均值向量, $\Sigma \in S_{++}^{n \times n}$ 半正定实对称矩阵为 $\vec{X}$ 的协方差矩阵, 则存在满秩矩阵 $B \in R^{n \times n},$ 使得 $\vec{Z}=B^{-1}(\vec{X}-\vec{\mu}),$ 而 $\vec{Z} \sim \mathcal{N}(\overrightarrow{0},\mathbf{I})$.**

（满秩的原因是我们需要求逆）

有了定理1, 我们就可以对随机向量$\vec{X}$ 做相应的线性变换, 使其随机变量在线性变换后彼此独立, 从而求出其联合概率密度函数, 具体地


$$
\begin{aligned}
\because \vec{Z}&=B^{-1}(\vec{X}-\vec{\mu}), \vec{Z} \sim \mathcal{N}(\overrightarrow{0}, I)\\
\therefore \quad p\left(z_{1}, \cdots, z_{n}\right) &=\frac{1}{(2 \pi)^{\frac{n}{2}}} \cdot e^{-\frac{1}{2} \cdot\left(Z^{\top} Z\right)} \\
p\left(z_{1}\left(x_{1}, \cdots, x_{n}\right), \cdots\right) &=\frac{1}{(2 \pi)^{\frac{n}{2}}} \cdot e^{-\frac{1}{2} \cdot\left[\left(B^{-1}(\vec{X}-\vec{\mu})\right)^{\top}\left(B^{-1}(\vec{X}-\vec{\mu})\right)\right]} \\
&=\frac{1}{(2 \pi)^{\frac{n}{2}}} \cdot e^{-\frac{1}{2} \cdot\left[(\vec{X}-\vec{\mu})^{\top}\left(B B^{\top}\right)^{-1}(\vec{X}-\vec{\mu})\right]} 
\end{aligned}
\tag{7}
$$
$$
\begin{aligned}
\therefore \quad 1 &=\int_{-\infty}^{+\infty} \cdots \int_{-\infty}^{+\infty} p\left(z_{1}\left(x_{1}, \cdots, x_{n}\right), \cdots\right) d z_{1} \cdots d z_{n} \\
&=\int_{-\infty}^{+\infty} \cdots \int_{-\infty}^{+\infty} \frac{1}{(2 \pi)^{\frac{n}{2}}} \cdot e^{-\frac{1}{2} \cdot\left[(\vec{X}-\vec{\mu})^{\top}\left(B B^{\top}\right)^{-1}(\vec{X}-\vec{\mu})\right]} d z_{1} \cdots d z_{n}
\end{aligned}
\tag{8}
$$



由多元函数换元变换公式, 我们还需要求出雅可比行列式 $J(\vec{Z} \rightarrow \vec{X}),$ 由(7)可得
$$
J(\vec{Z} \rightarrow \vec{X})=\left|B^{-1}\right|=|B|^{-1}=|B|^{-\frac{1}{2}} \cdot\left|B^{\top}\right|^{-\frac{1}{2}}=\left|B B^{\top}\right|^{-\frac{1}{2}}
\tag{9}
$$
由(8)(9), 我们可进一步得
$$
1=\int_{-\infty}^{+\infty} \cdots \int_{-\infty}^{+\infty} \frac{1}{(2 \pi)^{\frac{n}{2}}\left|B B^{\top}\right|^{\frac{1}{2}}} \cdot e^{-\frac{1}{2} \cdot\left[(\vec{X}-\vec{\mu})^{\top}\left(B B^{\top}\right)^{-1}(\vec{X}-\vec{\mu})\right]} d x_{1} \cdots d x_{n}
\tag{10}
$$
我们得到随机向量$\vec{X} \sim \mathcal{N}(\vec{\mu}, \Sigma)$ 的联合概率密度函数为
$$
p\left(x_{1}, \cdots, x_{n}\right)=\frac{1}{(2 \pi)^{\frac{n}{2}}\left|B B^{\top}\right|^{\frac{1}{2}}} \cdot e^{-\frac{1}{2} \cdot\left[(\vec{X}-\vec{\mu})^{\top}\left(B B^{\top}\right)^{-1}(\vec{X}-\vec{\mu})\right]}
\tag{11}
$$
在(11)中, 随机向量$\vec{X}$ 的协方差矩阵还未得到体现, 我们可通过线性变换(7)做进一步处理
$$
\begin{aligned}
\Sigma &=\mathbf{E}\left[(\vec{X}-\vec{\mu})(\vec{X}-\vec{\mu})^{\top}\right] \\
&=\mathbf{E}\left[(B \vec{Z}-\overrightarrow{0})(B \vec{Z}-\overrightarrow{0})^{\top}\right] \\
&=\operatorname{Cov}(B \vec{Z}, B \vec{Z}) \\
&=B \operatorname{Cov}(\vec{Z}, \vec{Z}) B^{\top} \\
&=B B^{\top}
\end{aligned}
\tag{12}
$$
>  常用结论：$Cov(B\vec{Z})=BCov(\vec{Z})B^{\top}$

我们发现, (11)中 $B B^{\top}$ 就是线性变换前的随机向量 $\vec{X} \sim \mathcal{N}(\vec{\mu}, \Sigma)$ 的协方差矩阵 $\Sigma$, 所以由(11)
(12), 我们可以得到联合概率密度函数的最终形式
$$
p\left(x_{1}, \cdots, x_{n}\right)=\frac{1}{(2 \pi)^{\frac{n}{2}}|\Sigma|^{\frac{1}{2}}} \cdot e^{.-\frac{1}{2} \cdot [(\vec{X}-\vec{\mu})^{\top} \Sigma^{-1}(\vec{X}-\vec{\mu})]}
\tag{13}
$$
这就是多元高斯分布的概率密度函数，非常重要，建议刻进DNA里。

