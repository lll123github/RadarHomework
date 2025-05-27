非常精彩的问题！你现在关心的是：

> **根据阵列因子 $AF(\theta, \phi)$ 的表达式，它的最大值与阵元数 $N_{\text{az}}, N_{\text{el}}$ 之间是什么关系？**

这涉及到阵列方向图主瓣的**归一化与峰值**，也是雷达系统中评估天线性能（如最大增益）时的核心参数。

---

## ✅ 一、阵列因子（AF）的理论最大值

我们回顾二维均匀矩形阵列（URA）中阵列因子的一般表达式（忽略相位扫描项）：

$$
AF(\theta, \phi) = \sum_{m=0}^{N_{\text{az}}-1} \sum_{n=0}^{N_{\text{el}}-1} e^{j k d \left[ m \sin\theta \cos\phi + n \sin\phi \right]}
$$

当信号从阵列的**正前方**入射（即主瓣方向：$\theta = 0, \phi = 0$）时：

* $\sin\theta = 0$
* $\sin\phi = 0$

此时：

$$
AF(0, 0) = \sum_{m=0}^{N_{\text{az}}-1} \sum_{n=0}^{N_{\text{el}}-1} 1 = N_{\text{az}} \cdot N_{\text{el}}
$$

---

## ✅ 二、解析表达式中如何体现这一点

你还可以从**解析表达式**角度看它：

$$
AF(\theta, \phi) =
\left[
\frac{\sin\left( N_{\text{az}} \frac{\pi d}{\lambda} \sin\theta \cos\phi \right)}{N_{\text{az}} \sin\left( \frac{\pi d}{\lambda} \sin\theta \cos\phi \right)}
\right]
\cdot
\left[
\frac{\sin\left( N_{\text{el}} \frac{\pi d}{\lambda} \sin\phi \right)}{N_{\text{el}} \sin\left( \frac{\pi d}{\lambda} \sin\phi \right)}
\right]
\cdot (N_{\text{az}} \cdot N_{\text{el}})
$$

所以你可以这样理解：

* 每一项是一个 **normalized Dirichlet sinc**，即归一化阵列方向图；
* **最大值出现在 $\theta = 0, \phi = 0$ 时**，此时两个正弦项都趋于 0，根据极限性质：

$$
\lim_{x \to 0} \frac{\sin(N x)}{N \sin(x)} = 1
$$

因此：

$$
AF_{\max} = N_{\text{az}} \cdot N_{\text{el}}
$$

---

## ✅ 三、物理意义与增益关系

* 这是阵列因子在主瓣方向（正前方）叠加的最大幅度；
* 如果你把这个幅度平方，得到的是功率方向图峰值：

  $$
  |AF_{\max}|^2 = (N_{\text{az}} \cdot N_{\text{el}})^2
  $$
* 这个值正比于**最大方向性增益（不含效率）**：

  $$
  G_D = |AF_{\max}|^2 \Big/ \int_{\Omega} |AF(\theta,\phi)|^2 d\Omega
  $$

---

## ✅ 四、总结

| 项目           | 关系                                              |            |                                          |
| ------------ | ----------------------------------------------- | ---------- | ---------------------------------------- |
| 峰值方向（主瓣中心）   | $\theta = 0, \phi = 0$                          |            |                                          |
| 幅度最大值        | $AF_{\max} = N_{\text{az}} \cdot N_{\text{el}}$ |            |                                          |
| 功率最大值（方向图峰值） | (                                               | AF\_{\max} | ^2 = (N\_{\text{az}} N\_{\text{el}})^2 ) |
| 与方向性增益关系     | 是方向性和增益的直接来源                                    |            |                                          |

---

如果你正在构造一个方向图函数 $G(\theta, \phi)$，其最大值也应匹配这个结果：

$$
G_{\max} \propto |AF_{\max}|^2 = (N_{\text{az}} N_{\text{el}})^2
$$

是否需要我进一步帮你将这个最大值归一化，构造真实的方向图 $G(\theta, \phi)$ 函数？
