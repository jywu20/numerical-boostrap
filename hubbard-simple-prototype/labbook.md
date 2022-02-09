复现2006.06002中对Hubbard模型的处理
======

Created on 2022.1.17.

我们这里尝试复现2006.06002中对Hubbard模型的处理，并且为更多的模型提供一个原型。

模型是
$$
H=-\sum_{\langle x y\rangle \sigma} c_{x \sigma}^{\dagger} c_{y \sigma}+U \sum_{x} n_{x \uparrow} n_{x \downarrow}
$$
在二维正方晶格上面。

# 理论

- Eq. (2)中的$C_\alpha$和$U_\alpha$包括：
  - 连续对称性：总粒子数、总自旋
  - 离散对称性：动量、旋转
- $\mathcal{C}_2$集合（保证正定性）
  - 由于局域性，由足够短的费米算符串张成
- $\mathcal{C}_1$集合（保证共轭、连续对称性生成元守恒、离散对称性成立）
  - 由$\mathcal{C}_2$中两个串相乘得到的全体算符张成
- 应当在Eq. (2)的约束下优化$F[H]$，这能够给系统总能量设置一个下界（是下界，因为真实的能量肯定在可行域中）

# 算法纲要

## $\mathcal{C}_2$的表示

为了避免double counting，应当始终使用正规序。

2022.1.17 各点的排序方式应该是怎样的？

## 泛函$\mathcal{F}$的表示

泛函$\mathcal{F}[O] := \langle \rho O \rangle$是一个线性泛函。按照前述$C_1$和$C_2$的构造，我们可以将$\mathcal{F}$看成线性映射$\mathcal{C}_2 \times \mathcal{C}_2 \to \Complex$。
最为朴素的实现是将$\mathcal{F}$当成一个列矢量，其基底是$\mathcal{C}_2 \times \mathcal{C}_2$。

Poor man's OPE: 用反对易关系把一个算符串重排成一系列正规序算符的线性组合然后代入$\mathcal{F}$中；一个费米算符出现两次就认为是零（因为不可能有两个一样的费米子）。

## 晶格

由于在无穷大的晶格上工作（远处的算符平移到原点附近来求值），不使用任何基于有限晶格的库。看起来，将每个格点都标上号是没有意义的，于是我们使用坐标$(x, y)$来标记格点。

哈密顿量中的各项也是平移到原点附近求值，最后整个哈密顿量正比于格点总数$N$，把$N$扔掉以后对能量做最优化。
$$
H=-\sum_{\langle x y\rangle \sigma} c_{x \sigma}^{\dagger} c_{y \sigma}+U \sum_{x} n_{x \uparrow} n_{x \downarrow}
$$
的第一项正比于$2N$，第二项正比于$N$，即
$$
\begin{aligned}
    \langle H \rangle &= - 2 N t \sum_\sigma \langle c^\dagger_{0 \sigma} c_{\hat{\boldsymbol{x}} \sigma} \rangle + U N \lang n_{0\uparrow} n_{0\downarrow} \rang  \\
    &= - 2 N t \sum_\sigma \langle c^\dagger_{0 \sigma} c_{\hat{\boldsymbol{x}} \sigma} \rangle + U N \lang c^\dagger_{0 \uparrow} c^\dagger_{0 \downarrow} c_{0 \downarrow} c_{0 \uparrow} \rang .
\end{aligned}
$$
类似的，
$$
\lang N \rangle = N \sum_\sigma \lang c^\dag_{0 \sigma} c_{0 \sigma} \rang,
$$
以及
$$
\lang S_z \rang = N (\lang c^\dag_{0 \uparrow} c_{0 \uparrow} \rang - \lang c^\dag_{0 \downarrow} c_{0 \downarrow} \rang).
$$

## 施加对称性

有待解决的问题：是否有把握
- 或者所有约束条件中涉及的算符序列都可以通过Poor man's OPE划归到$\mathcal{C}_2$中，或是
- 所有约束条件都能够化归为对$\mathcal{F}$的约束？看起来和$\mathcal{C}_1$有关的约束主要是用于构造$\mathcal{F}$的；主要的计算量都在正定性条件上面。

- 平移对称性：
- 旋转对称性：
- 复共轭：似乎可以将$\mathcal{C}_2$集合中互为共轭转置的算符配对

## 正定性

我怀疑这才是真的关键的地方……这里的约束看起来是不能够通过限制$\mathcal{F}$的形式引入的

2022.1.17 感觉应该先做一个toy model感觉一下情况；应该看一下2004.10212那篇文章

2022.1.30 目前的计划：
- 使用$\lang c^\dag_{i \sigma} c_{j \sigma} \rang$和$\lang n_{i \uparrow} n_{i \downarrow} \rang$做变量
- 优化目标$E$的表达式是显然的
- 通过$\lang H \mathcal{O} \rang$可以制造出适当的OPE
- 这样的问题是，无法在$\lang c^\dag_{i \sigma} c_{j \sigma} \rang$和$\lang n_{i \uparrow} n_{i \downarrow} \rang$之间建立联系……
  - 我需要一个Hubbard模型版本的位力定理，
  - 或者说要仅仅使用$\lang c^\dag_{i \sigma} c_{j \sigma} \rang$和$\lang n_{i \uparrow} n_{i \downarrow} \rang$构造一个方程

最好是拿mathematica算一算对易子

# 实验记录

## 2022.2.8

在`hubbard-simple-prototype\2022-2-8.nb`中计算费米子对易子。

下面要写一个程序弄一个类似于`hubbard-simple-prototype\2d-square-lattice-indexing.pdf`的格点编码。

## 2022.2.9

完成格点编码。准备构造$\mathcal{O}_2$集合。在`hubbard-simple-prototype\2022-2-9.nb`中做这件事。