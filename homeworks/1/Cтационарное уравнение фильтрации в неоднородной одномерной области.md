Для внутренних ячек:
$$
\begin {equation}
y_i\Delta x \sim f(x_{i + \frac{1}{2}}) \cdot \Delta x = q_{i+1} - q_i 
= - \frac{2K_{i+1}K_{i}}{K_i+K_{i+1}}\cdot\frac{h_{i+1} - h_i}{\Delta x}
+ \frac{2K_{i-1}K_{i}}{K_i+K_{i-1}}\cdot\frac{h_i - h_{i-1}}{\Delta x}=\\

\end {equation}
$$
$$
\begin {equation}
= \frac{2K_i}{\Delta x}
(
- h_{i-1}\frac{K_{i-1}}{K_i+K_{i-1}}
- h_{i+1}\frac{K_{i+1}}{K_i+K_{i+1}}
+ h_i\frac{K_{i+1}(K_i+K_{i-1}) + K_{i-1}(K_i+K_{i+1})}{(K_i+K_{i+1})(K_i+K_{i-1})}
)
\end {equation}
$$
$$
\begin {equation}
\Rightarrow 

- h_{i-1}\frac{K_{i-1}}{K_i+K_{i-1}}
- h_{i+1}\frac{K_{i+1}}{K_i+K_{i+1}}
+ h_i\frac{K_{i+1}(K_i+K_{i-1}) + K_{i-1}(K_i+K_{i+1})}{(K_i+K_{i+1})(K_i+K_{i-1})}
= f(x_{i + \frac{1}{2}}) \cdot \frac{(\Delta x)^2}{2K_i}
\end {equation}
$$

Для левой ячейки:
$$
\begin {equation}
y_0\Delta x \sim f(x_{\frac{1}{2}}) \cdot \Delta x 
= - \frac{2K_{1}K_{0}}{K_0+K_{1}}\cdot\frac{h_{1} - h_0}{\Delta x}
+ \frac{2K_{-\frac{1}{2}}K_{0}}{K_0+K_{-\frac{1}{2}}}\cdot\frac{h_0 - a}{\Delta x/2}=\\

\end {equation}
$$
$$
\begin {equation}
\Rightarrow
- \frac{K_{1}}{K_0+K_{1}}\cdot h_{1}
+ \frac{2K_{-\frac{1}{2}}(K_0 + K_1)+K_1(K_0+K_{-\frac{1}{2}})}{(K_0+K_{-\frac{1}{2}})(K_0+K_1)}\cdot h_0
= \frac{(\Delta x)^2}{2K_0} \cdot f(x_{\frac{1}{2}}) + \frac{2a}{K_0+K_{-\frac{1}{2}}}
\end {equation}
$$
Аналогично для правой ячейки:

$$
\begin {equation}
y_{n-1}\Delta x \sim f(x_{n-\frac{1}{2}}) \cdot \Delta x 
= - \frac
    {2K_{n-1}K_{n-\frac{1}{2}}}{K_{n-1}+K_{n-\frac{1}{2}}}\cdot\frac{b - h_{n-1}}
    {\Delta x/2}
+ \frac{2K_{n-1}K_{n-2}}{K_{n-1}+K_{n-2}}\cdot\frac{h_{n-1} - h_{n-2}}{\Delta x}=\\
\end {equation}
$$
$$
\begin {equation}
  \frac
    {2K_{n-\frac{1}{2}}(K_{n-1}+K_{n-2})+K_{n-2}(K_{n-1}+K_{n-\frac{1}{2}})}
    {(K_{n-1}+K_{n-\frac{1}{2}})(K_{n-1}+K_{n-2})} \cdot h_{n-1}
- \frac{K_{n-2}}{K_{n-1}+K_{n-2}}\cdot h_{n-2}
= \frac{(\Delta x)^2}{2K_{n-1}} \cdot f(x_{n-\frac{1}{2}}) 
+ \frac{2b}{K_{n-1}+K_{n-\frac{1}{2}}}
\end {equation}
$$