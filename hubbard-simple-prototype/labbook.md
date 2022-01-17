复现2006.06002中对Hubbard模型的处理
======

我们这里尝试复现2006.06002中对Hubbard模型的处理，并且为更多的模型提供一个原型。

$$
H=-\sum_{\langle x y\rangle \sigma} c_{x \sigma}^{\dagger} c_{y \sigma}+U \sum_{x} n_{x \uparrow} n_{x \downarrow}
$$
在正方晶格上面。

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

## 泛函$\mathcal{F}$的表示

泛函$\mathcal{F}[O] := \langle \rho O \rangle$是一个线性泛函。按照前述$C_1$和$C_2$的构造，我们可以将$\mathcal{F}$看成线性映射$\mathcal{C}_2 \times \mathcal{C}_2 \to \Complex$。
最为朴素的实现是将$\mathcal{F}$当成一个列矢量，其基底是$\mathcal{C}_2 \times \mathcal{C}_2$。

## 晶格



# 实验记录
