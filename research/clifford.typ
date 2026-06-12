#import "@preview/physica:0.9.7":*
#set text(font:"Noto Sans CJK SC", size:12pt)
*例子*

考虑波函数
$
  ket(psi) = 1/sqrt(2) * (ket(0)^(times.circle N) + ket(1)^(times.circle N))
$

这个函数很难被local update sample！

如果进行Clifford 演化呢?

$
  ket("Stabilizer"_i) = cal(C) ket(0)^(times.circle N)
$
其中 C 是一个Clifford 门阵列。

这样采样的波函数就变成了所有被Clifford Group stabilize的集合

直觉上：这个态比原来的GHZ态好采样很多！

进一步地，我们可以证明，Clifford演化的终态是最适合采样的（最均匀分布在计算基上）

== 逆参与率 (IPR)
对于一个在计算基 ${ket(x)}$ 下展开的量子态 $ket(psi) = sum_(x=1)^D psi_x ket(x)$，其逆参与率定义为：
$ "IPR"(ket(psi)) = sum_(x=1)^D abs(psi_x)^4 $
其中归一化条件为 $sum abs(psi_x)^2 = 1$。

== 极值条件
- *最大值*：当态完全局域化（如 $ket(x)$）时，$"IPR" = 1$。
- *最小值*：当态完全离域（均匀叠加，Flat state）时，即 $abs(psi_x) = 1/sqrt(D)$，
  $ "IPR"_min = sum_(x=1)^D (1/sqrt(D))^4 = D dot 1/D^2 = 1/D $

*目标*：证明随机 Clifford 演化产生的态，其 IPR 接近 $1/D$。

= 证明过程

证明分为两个视角：*结构视角（单次演化）*和*统计视角（平均演化）*。
#let tensor = $times.circle$
== 结构视角：Clifford 态的“平坦性”

Clifford 群 $cal(C)_n$ 作用于基态 $ket(0)^(tensor n)$ 产生的态被称为*稳定子态 (Stabilizer States)*。

*引理 1：稳定子态是平坦的 (Flat)。*

任何稳定子态 $ket(psi)$ 在计算基下的振幅分布只有两种可能：要么为 0，要么模长相等。即：
$ abs(psi_x) in { 0, 1/sqrt(K) } $
其中 $K$ 是波函数非零分量（Support）的大小，$1 <= K <= D$。

*推论：*

对于一个稳定子态，其 IPR 严格等于其支撑集大小的倒数：
$ "IPR" = sum_(x: psi_x != 0) (1/sqrt(K))^4 = K dot 1/K^2 = 1/K $

这意味着，*只要 Clifford 演化使波函数扩散到整个希尔伯特空间（即 $K -> D$），它就能达到数学上的绝对最小值 $1/D$*。与高斯波包或其他非平坦分布不同，Clifford 态没有“拖尾”，这使得它在同等扩散程度下 IPR 最小。

== 统计视角：酉 2-Design 性质

为了证明随机演化确实会使 $K$ 变大（即 IPR 变小），我们利用 Clifford 群的统计性质。

*定理：Clifford 群构成一个酉 2-design (Unitary 2-design)。*

这意味着对于任何算符多项式达二阶矩的平均值，Clifford 群上的平均等同于 Haar 测度（完全随机幺正群）上的平均。

数学表述为：
$ bb(E)_(U in "Clifford") [ (U ketbra(psi_0, psi_0) U^dagger)^(tensor 2) ] = bb(E)_(U in "Haar") [ (U ketbra(psi_0, psi_0) U^dagger)^(tensor 2) ] $

*计算平均 IPR：*

我们计算经过随机 Clifford 演化后态 $ket(psi) = U ket(0)$ 的平均 IPR。由于 IPR 是波函数振幅的四次函数（即密度矩阵的二次函数 $tr(rho^2)$ 的某种形式，具体为 $sum abs(psi_x)^4$），它仅依赖于二阶矩。

利用 2-design 性质，我们可以直接使用 Haar 随机态的已知结果。对于希尔伯特空间维数为 $D$ 的 Haar 随机态，平均 IPR 为：
$ bb(E)["IPR"] = 2/(D+1) $

*推导简述：*

对于随机态，$psi_x$ 服从 Porter-Thomas 分布。
$ bb(E)[abs(psi_x)^4] = 2 (bb(E)[abs(psi_x)^2])^2 = 2 (1/D)^2 = 2/D^2 $
求和所有 $D$ 个分量：
$ bb(E)["IPR"] = sum_(x=1)^D bb(E)[abs(psi_x)^4] = D dot 2/D^2 = 2/D $
#text(size: 0.8em, style: "italic")[(注：精确公式为 $2/(D+1)$，但在大 $D$ 极限下 $approx 2/D$)]。

= 结论

+ *理论最小值*：IPR 的绝对下界是 $1/D$。
+ *Clifford 结果*：随机 Clifford 演化产生的态，其 IPR 的期望值为 $2/(D+1) approx 2/D$。
+ *对比*：$2/D$ 与最小值 $1/D$ 同阶，且随着系统尺寸指数衰减。

*综上所述，Clifford 随机演化通过其 2-design 的混沌特性 (Scrambling)，能够将初态演化为希尔伯特空间中的典型态。这些态在结构上是平坦的，在统计上其 IPR 达到了物理允许的最小值量级，证明了其最大化离域的能力。*