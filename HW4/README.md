# Homework 4: Selection-Model Marginal Effects and Panel Data

📄 [Problem set](assignment/HW4-questions.pdf) · 📘 [Full solutions (PDF)](report/HW4-solutions.pdf) · 💻 [Code](code/)

| Question | Topic | Script |
|:--:|---|---|
| Q1 | Marginal effects in the sample-selection model | analytical (see report) |
| Q2 | Fixed-effects (within) estimator: assumptions and derivation | analytical (see report) |
| Q3 | Do executions deter murder? Pooled OLS, first differences, FE | [`Q3_murder_panel.m`](code/Q3_murder_panel.m) |

Run with `run('HW4/code/run_all.m')`.

---

## Q1: Marginal effects in the selection model

With $E[y_2\mid y_1^*>0,x]=x_2'\beta_2+\sigma_{12}\lambda(\omega_1)$, where $\omega_1=x_1'\beta_1$, $\lambda=\phi/\Phi$, and $\lambda'(\omega)=-\lambda(\omega)\,[\omega+\lambda(\omega)]$:

| Covariate appears in | Marginal effect at $\bar x$ |
|---|---|
| participation equation only | $-\beta_{1,j}\,\sigma_{12}\,\lambda(\bar\omega_1)\big[\bar\omega_1+\lambda(\bar\omega_1)\big]$ |
| outcome equation only | $\beta_{2,j}$ |
| both equations | $\beta_{2,j}-\beta_{1,j}\,\sigma_{12}\,\lambda(\bar\omega_1)\big[\bar\omega_1+\lambda(\bar\omega_1)\big]$ |

## Q2: The within estimator

Demeaning $y_{it}=x_{it}'\beta+c_i+\varepsilon_{it}$ with $M_T=I_T-\tfrac1T\mathbf 1_T\mathbf 1_T'$ removes $c_i$:

$$\tilde y=(I_n\otimes M_T)\,y=\tilde X\beta+\tilde\varepsilon,\qquad \hat\beta_{FE}=(\tilde X'\tilde X)^{-1}\tilde X'\tilde y=\Big(\sum_{i,t}\tilde x_{it}\tilde x_{it}'\Big)^{-1}\sum_{i,t}\tilde x_{it}\tilde y_{it}$$

$$\hat c_i=\bar y_i-\bar x_i'\hat\beta_{FE},\qquad \hat\sigma^2_\varepsilon=\frac{1}{n(T-1)-k}\sum_{i,t}\hat{\tilde\varepsilon}_{it}^{\,2}$$

## Q3: Executions and murder rates (US states)

$$\text{mrdrte}_{it}=\eta_t+\beta_1\,\text{exec}_{it}+\beta_2\,\text{unemp}_{it}+c_i+\varepsilon_{it}$$

If executions deter murder, we expect β₁ < 0. Unemployment is expected to raise crime, so β₂ > 0.

| Specification | β̂ exec | SE | p-value |
|---|---:|---:|---:|
| Pooled OLS, 1990 & 1993 (HC1) | 0.1277 | 0.1342 | 0.344 |
| First difference, 1990 → 1993 | **−0.0929** | 0.0160 | < 0.001 |
| First difference, **excluding Texas** | −0.0549 | 0.0770 | 0.480 |
| Within (FE), 1987/90/93, SEs clustered by state | −0.1383 | 0.0787 | 0.079 |

**Takeaways**

- Pooled OLS shows no deterrent effect. It ignores the state effects $c_i$, which are correlated with execution policy.
- Removing $c_i$ by first differencing gives a negative and highly significant coefficient.
- **Texas drives that result.** Texas had 34 executions in 1991–93, 23 more than the next state, Virginia, with 11. Without Texas the effect shrinks and loses significance.
- With all three years and cluster-robust inference, the effect is negative but only marginally significant. The evidence for deterrence is fragile.
