
# Varied-cost markets

``
\begin{align}
\text{maximize}\quad & v(\mathcal{X}) =  \operatorname{E}\Bigl[\max\bigr\{t_0,
\max\{t_j Z_j : j \in \mathcal{X}\}\bigr\}\Bigr] \\
\text{subject to}\quad & \mathcal{X} \subseteq \mathcal{C}, ~~\sum_{j\in \mathcal{X}} g_j \leq H
\end{align}
``


```@docs
VariedCostsMarket
```

```@ref
valuation
```
