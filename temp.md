图中显示的函数为：

$$
\frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt
$$

这是误差函数（Error Function）的定义，其在 MATLAB 中由内建函数 `erf(x)` 表示。

---

### MATLAB 中的转换方式

在 MATLAB 中，您无需自己计算该积分，可以直接使用：

```matlab
y = erf(x);
```

这个函数 `erf(x)` 的定义就是：

$$
\text{erf}(x) = \frac{1}{\sqrt{2\pi S/N}} \int_{\beta}^{\inf} e^{\frac{-y^2}{2 S/N}} dy
$$

因此，**图中所示函数即为 MATLAB 的 `erf(x)` 函数本身**，不需要进一步转换。只要你在 MATLAB 中使用 `erf(x)`，调用的正是图中表达的数学函数。

---

如果您还需要与其他函数（如 `erfc(x)`，互补误差函数）进行比较或结合使用，我也可以提供帮助。是否需要展示这部分？
