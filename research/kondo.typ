#import "@preview/physica:0.9.7": *

#set math.equation(numbering: "(1)")

= Yao-Lee-Kondo 模型的 Majorana 重新推导

#let up = $arrow.t$
#let down = $arrow.b$
#let vec(body) = $arrow(body)$

== 1. 传导电子算符定义
将传导电子分成四个majorana算符
$
  c_(j up) = 1/2 (eta^j_(1 up) - i eta^j_(2 up)), quad c_(j down) = 1/2 (eta^j_(1 down) - i eta^j_(2 down)) \
  c_(j up)^dagger = 1/2 (eta^j_(1 up) + i eta^j_(2 up)), quad c_(j down)^dagger = 1/2 (eta^j_(1 down) + i eta^j_(2 down))
$ <eq-convention>

在此约定下，占据数算符为：
$ n_(j up) = c_(j up)^dagger c_(j up) = 1/4 (eta^j_(1 up) + i eta^j_(2 up))(eta^j_(1 up) - i eta^j_(2 up)) = 1/2 - i/2 eta^j_(1 up) eta^j_(2 up) $

== 2. 传导电子自旋 $vec(s)_j$ 的推导

利用自旋算符标准定义 $s^a = 1/2 c^dagger sigma^a c$：

- *Z 分量*:
$ s_j^z &= 1/2 (n_(j up) - n_(j down)) \
      &= 1/2 [(1/2 - i/2 eta^j_(1 up) eta^j_(2 up)) - (1/2 - i/2 eta^j_(1 down) eta^j_(2 down))] \
      &= - i / 4 (eta^j_(1 up) eta^j_(2 up) - eta^j_(1 down) eta^j_(2 down)) $

- *X 分量*:
$ s_j^x &= 1/2(S^+_j + S^-_j ) \
      &=1/2 (c_(j up)^dagger c_(j down) + c_(j down)^dagger c_(j up)) \
      &= 1/2 [ 1/4 (eta^j_(1 up) + i eta^j_(2 up))(eta^j_(1 down) - i eta^j_(2 down)) + 1/4 (eta^j_(1 down) + i eta^j_(2 down))(eta^j_(1 up) - i eta^j_(2 up)) ] \
      &= 1/8[{eta^j_(1 up), eta^j_(1 down)} + {eta^j_(2 up), eta^j_(2 down)} + i [eta^j_(2 up), eta^j_(1 down)] - i [eta^j_(1 up), eta^j_(2 down)]] \
      &= - i / 4 (eta^j_(1 up) eta^j_(2 down) - eta^j_(2 up) eta^j_(1 down)) $

- *Y 分量*:
$ s_j^y &= 1/(2i) (S^+_j - S^-_j ) \
      &=1/(2i) (c_(j up)^dagger c_(j down) - c_(j down)^dagger c_(j up)) \
      &= 1/(8i) [(eta^j_(1 up) + i eta^j_(2 up))(eta^j_(1 down) - i eta^j_(2 down)) - (eta^j_(1 down) + i eta^j_(2 down))(eta^j_(1 up) - i eta^j_(2 up))] \
      & = 1/(8i) [[eta^j_(1 up), eta^j_(1 down)] - i {eta^j_(1 up), eta^j_(2 down)} + i {eta^j_(2 up), eta^j_(1 down)} +[eta^j_(2 up),eta^j_(2 down)]] \

      &= - i / 4 (eta^j_(1 up) eta^j_(1 down) + eta^j_(2 up) eta^j_(2 down)) $

== 3. Yao-Lee 局域自旋 $vec(S)_j$ 的表示

为了保持 $S=1/2$ 代数且满足 $chi^2=1$，局域自旋定义为：
$ arrow(S) = -i/2 arrow(chi) times arrow(chi) $

$ S_j^x = -i / 2 chi_j^y chi_j^z, quad S_j^y = -i / 2 chi_j^z chi_j^x, quad S_j^z = -i / 2 chi_j^x chi_j^y $

== 4. 最终近藤耦合项 $H_K$

将 $vec(s)_j$ 与 $vec(S)_j$ 代入 $H_K = J sum_j vec(S)_j dot vec(s)_j$。由于 $S^a$ 和 $s^a$ 均带有 $-i$ 系数，其乘积项出现 $(-i)^2 = -1$：

$ H_K = - J / 8 sum_j [   &chi_j^x chi_j^y (eta^j_(1 up) eta^j_(2 up) - eta^j_(1 down) eta^j_(2 down)) \
  &+ chi_j^y chi_j^z (eta^j_(1 up) eta^j_(2 down) - eta^j_(2 up) eta^j_(1 down)) \
  &+ chi_j^z chi_j^x (eta^j_(1 up) eta^j_(1 down) + eta^j_(2 up) eta^j_(2 down)) ] $



// == HS transform for Majorana
// Take one term out,
// let
// $ H_(K 1) =J/8 sum_j (i chi_j^x chi_j^y) (i eta^j_(1 up) eta^j_(2 up)) $

// $ exp(- Delta tau H_(K 1)) = (1 / 2)^L sum_( {s_j} = plus.minus 1) exp[ sum_j ( lambda s_j (i chi_j^x chi_j^y + i eta^j_(1 up) eta^j_(2 up)) + (Delta tau J) / 8 ) ] $ 

// = Sign problem



