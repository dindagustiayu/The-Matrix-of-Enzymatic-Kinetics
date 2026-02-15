[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://github.com/dindagustiayu/The-Matrix-of-Enzymatic-Kinetics/blob/main/Michaelis-Menten%20model.cpp)

# The Michael-Menten model and Linear-Burk Plot

## Michael-Menten Kinetics
Michaelis-Menten kinetics describes the rate of enzymatic reactions, providing a model to understand how enzymes interact with substrates. The general reaction scheme of an enzyme-catalyzed reaction is as follows:
<p align='center'>
    $$E + S\; \overset{k_1}{\underset{k_{-1}}{\longrightleftharpoons}} \;ES\; \overset{k_2}{\longrightarrow}\;P + E$$
</p>
The enzyme interacts with the substrate by binding to its active site to form the enzyme-substrate complex, ES. That reastion is followed by the decomposition of ES to generate the free enzyme, E, and the new product, P.

| Rate constant | Reaction |
| ------------- | -------- |
| $k_{1}$ | The binding of the enzyme to the substrate, forming the enzyme-substrate complex. |
| $k_{2}$ | Catalytic rate; the catalysis reaction producing the final reaction product and regenerating the free enzyme. This is the rate-limiting step. |

## Lineweaver-Burk Plot
For many such catalyzed reactions, it is found that:
<p align='center'>
    $$v_{0}=\frac{(V_{max}[S])}{(K_{M}+[S])}$$
    $$\frac{1}{v}=\frac{(K_M + [S])}{v_{max}[S]}$$
    $$\frac{1}{v}=\left(\frac{K_M}{v_{max}} \right)\;\left(\frac{1}{[S]} \right) + \frac{1}{v_{max}}$$
</p>
Apply this to equation for a straight line $y = mx + b$ and we have:
<p align='center'>
    $$y=\frac{1}{v}$$
    $$x=\frac{1}{[S]}$$
    $$m=slope=\frac{K_M}{V_{max}}$$
    $$b=y-intercept=\frac{1}{V_{max}}$$
</p>
When we plot $y=\frac{1}{v}$ versus $x=\frac{1}{[S]}$, we obtain a straight line.
<p align='center'>
    $$x-intercept=\frac{-1}{K_M}$$
    $$y-intercept=\frac{1}{V_{max}}$$
    $$slope=\frac{K_M}{V_{max}}$$
</p>

<p align='center'>
    <img src='https://upload.wikimedia.org/wikipedia/commons/thumb/4/49/LB-plot.jpg/500px-LB-plot.jpg' alt='Figure from Wekepedia'>
</p>


## P19.5 Exercise 
Pepsin is one of the  main enzymes in the digestive systems of mammals, where it helps to break down proteins into smallers peptides which can be absorbed by the small intestine. In a study of the enzymatic kinetics of pepsin with the protein bovine serum albumin (S), the rate of reaction, $v$, was measured as a function of substrate concentration, [S]; these data are given in the file [pepsin.txt](https://github.com/dindagustiayu/The-Matrix-of-Enzymatic-Kinetics/blob/main/pepsin.txt)

```
    # Reaction rates, v (in mM.s-1), for the enzymatic kinetics of pepsin with
    # bovine serum albumin (S, in mM)

    S / mM	v / mM.s-1
    0.1		0.00339
    0.2		0.00549
    0.5		0.00875
    1.0		0.01090
    5.0		0.01358
    10.0	0.01401
    20.0	0.01423
```

In the case of a straight-line regression ("simple linear regression") on n points using the model function $y=a+b_x$ as follow:
<p align='left'>
    $$\hat{\beta}=\begin{pmatrix} a\\b\end{pmatrix}\;where\;a=\frac{S_y S_{xx} - S_{xy} S_x}{nS_{xx} - S^2_x}\;and\;b=\frac{nS_{xy} S_y S_x}{nS_{xx} - S^2_x},$$
</p>
and the summary statistics are defined as:
<p align='center'>
    $$S_x=\overset{n}{\underset{i=1}{\sum x_i}},\;S_y=\overset{n}{\underset{i=1}{\sum Y_i}},\;S_{xx}=\overset{n}{\underset{i=1}{\sum x_i^2}},\;S_{xy}=\overset{n}{\underset{i=1}{\sum x_iy_i}}.$$
</p>

### Step 1. Plotted as a function of reaction rates for the enzymatic kinetics of pepsin
```python
import numpy as np
import matplotlib.pyplot as plt

# Open file data
S, v = np.genfromtxt('pepsin.txt', unpack=True, skip_header=3)

# Transform data
x = 1 / S   # 1/[S]
y = 1 / v   # 1/v

plt.figure(figsize=(6, 4))
plt.plot(x, y, '+', markersize=10)
plt.xlabel('$\\mathrm{1/[S]} \\; (mM)^{-1}$')
plt.ylabel('$\\mathrm{1/[v]} \\; s(mM)^{-1}$')
```
<p align='center'>
    <img src='https://github.com/dindagustiayu/The-Matrix-of-Enzymatic-Kinetics/blob/main/Linear%20fit%20svg/Plotted%20as%20a%20function%20of%20enzymatic%20kinetics.svg'>
</p>


### Step 2. Fitting a Line of Best fit

The experiment was carried out at $$35^{0}C$$ and a pH of 2 with a total concentration of pepsin of $$E_0 = 0.028 mM$$. Use `np.linalg.lstsq` to fit the data and obtain values for $$K_M, v_{max} , and K_2$$.

For simple linear regression on $n$ points a straight line model, $y=a+bx$, 

### Step 2. Fitting a Line of Best fit
The experiment were carried out at $35 ^{o}C$ and a pH of 2 with a totl concentration of pepsin of $[E]_{0}=0.028\;mM$. Use `np.linalg.lstsq` to fit the data and obtain values for $K_M,\;v_{max}\;and\;k_2$.

For simple linear regression on $n$ points a straight line model, $y=a+bx$, 
<p align='center'>
    $$\hat{\beta}=\begin{pmatrix} a\\ b \end{pmatrix}=(X^T X)^{-1} X^Ty$$
</p>

defines the best fit coefficients $a$ and $b$ where the matrices $y$ and $X$
<p align='center'>
    $$y=\begin{pmatrix} y_1\\y_2\\ \vdots \\y_n\end{pmatrix},\; and\; X=\begin{pmatrix}1 & x_1\\1 & x_2\\ \vdots \\1 & x_n\end{pmatrix}$$
</p>
We have
<p align='center'>
    $$X^T X=\begin{pmatrix} 1 & 1 & 1 & \cdots & 1\\ x_1 & x_2 & x_3 \cdots & x_n\end{pmatrix} \; \begin{pmatrix} 1 & x_1\\1 & x_2\\1 & x_3\\ \vdots & \vdots\\ 1 & x_n\end{pmatrix} = \begin{pmatrix} n & \sum{x_i}\\ \sum{x_i} & \sum{x_i^2}\end{pmatrix}$$
</p>
Taking the inverse,
<p align='center'>
    $$(X^TX X)^{-1}= \frac{1}{\Delta} \;\begin{pmatrix} \sum{x_i^2} & -\sum{x_i}\\ -\sum{x_i} & n\end{pmatrix}$$
</p>
where the determinant,
<p align='center'>
    $$\Delta=n\sum{x_i^2}-(\sum{x_i})^2 = nS_{xx} - S_{x}^2$$
</p>
Next,

<p align='center'>
    $$\begin{align*}(X^TX)^{-1}X^T&=\frac{1}{\Delta} \begin{pmatrix} \sum{x_i^2} & -\sum{x_i}\\ -\sum{x_i} & n\end{pmatrix} \begin{pmatrix} 1 & 1 & 1 & \cdots & 1\\ x_1 & x_2 & x_3 & \cdots & x_n\end{pmatrix} \\&=\frac{1}{\Delta} \begin{pmatrix} \sum{x_i^2}-x_1 \sum{x_1} & \sum{x_i^2}-x_2 \sum{x_i} & \cdots & \sum{x_i^2}-x_n \sum{xi}\\ -\sum{x_i}+ nx_1 & -\sum{x_i}+nx_2 & \cdots & -\sum{x_i}+nx_n \end{pmatrix} \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_n\end{pmatrix} \end{align*}$$
</p>

Finally,
<p align='center'>
    $$ \begin{align*}\hat{\beta}&= (X^T X)^{-1}X^Ty \\ &= \frac{1}{\Delta} \begin{pmatrix} \sum{x_i^2}-x_1 \sum{x_1} & \sum{x_i^2}-x_2 \sum{x_i} & \cdots & \sum{x_i^2}-x_n \sum{xi}\\ -\sum{x_i}+ nx_1 & -\sum{x_i}+nx_2 & \cdots & -\sum{x_i}+nx_n \end{pmatrix} \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_n\end{pmatrix}=\begin{pmatrix} a\\b \end{pmatrix} \end{align*}$$
</p>
where,
<p align='center'>
    $$\begin{align} a&= \frac{1}{\Delta} \left [y_1 \left(\sum{x_i^2}- x_1 \sum{x_i} \right)+y_2 \left(\sum{x_i^2}-x_2 \sum{x_i} \right)+ \cdots + y_n \left(\sum{x_i^2}-x_n\sum{x_i} \right) \right] \\&=\frac{1}{\Delta} \left[\left(\sum{y_i} \right) \left(\sum{x_i^2} \right)- \left(\sum{x_i y_i}\right) \left(\sum{x_i}\right) \right] \\ &=\frac{S_yS_{xx}-S_{xy}S_x}{nS_{xx}-S_x^2} \end{align}$$
</p>
and
<p align='center'>
    $$\begin{align} b&= \frac{1}{\Delta} \left[y_1 \left(nx_1-\sum{x_i} \right)+y_2 \left(nx_2-\sum{x_i} \right)+ \cdots + y_n \left(nx_n - \sum{x_i} \right) \right] \\&=\frac{1}{\Delta} \left [n \left(\sum{x_iy_i} \right)- \left(\sum{y_i} \right) \left(\sum{x_i} \right) \right]\\ &=\frac{nS_{xy}-S_x S_y}{nS_{xx}-S_x^2}\end{align}$$
    $$\begin{align} a&=\frac{1}{\Delta}[y_1(\sum{x_i^2}-x_1 \sum{x_i})+y_2(\sum{x_i^2}-x_2\sum{x_i})+ \cdots + y_n (\sum{x_i^2}-x_n\sum{x_i})] \\&=\frac{1}{\Delta}[(\sum{y_i}(\sum{x_i^2})-(\sum{x_i y_i})(\sum{x_i})] \\ &=\frac{S_yS_{xx}-S_{xy}S_x}{nS_{xx}-S_x^2} \end{align}$$
</p>
and
<p align='center'>
    $$\begin{align} b&= \frac{1}{\Delta}[y_1(nx_1-\sum{x_i})+y_2(nx_2-\sum{x_i})+ \cdots + y_n(nx_n - \sum{x_i})] \\&=\frac{1}{\Delta}[n(\sum{x_iy_i}-(\sum{y_i})(\sum{x_i})]\\ &=\frac{nS_{xy}-S_x S_y}{nS_{xx}-S_x^2}\end{align}$$
</p>

```python
# Design matrix
A = np.vstack([x, np.ones_like(x)]). T

# Least squares fit
m, b = np.linalg.lstsq(A, y, rcond=None)[0]

# Parameters Lineweaver-Burk Plot
vmax = 1 /b
kM = m * vmax
E0 = 0.028
k2 = vmax /E0

print(f"vmax = {vmax:.4f} mM/s")
print(f"Km   = {Km:.4f} mM")
print(f"k2   = {k2:.4f} s^-1")
```
```
vmax = 0.0145 mM/s
Km   = 0.3267 mM
k2   = 0.5166 s^-1
```
```python
#Plot Lineweaver-Burk Plot
x_fit = np.linspace(min(x), max(x), 100)
y_fit = m * x_fit + b

plt.figure(figsize=(6, 4))
plt.plot(x, y, '+', color='blue', markersize=8, label='Data')
plt.plot(x_fit, y_fit, color='black', label='Fit')
plt.xlabel('$\\mathrm{1/[S]} \\; (mM)^{-1}$')
plt.ylabel('$\\mathrm{1/[v]} \\; s(mM)^{-1}$')
plt.title('Lineweaver-Burk Plot for Pepsin')
plt.legend()
plt.show()
```
<p align='center'>
    <img src='https://github.com/dindagustiayu/The-Matrix-of-Enzymatic-Kinetics/blob/main/Linear%20fit%20svg/Linearweaver-Burk%20Plot%20for%20pepsin.svg'>
</p>

