# Homework 2: Multinomial Logit and Count Data

📄 [Problem set](assignment/HW2-questions.pdf) · 📘 [Full solutions (PDF)](report/HW2-solutions.pdf) · 💻 [Code](code/)

| Question | Topic | Script |
|:--:|---|---|
| Q1 | MNL derivation from Gumbel errors, own- and cross-marginal effects | analytical (see report) |
| Q2 | Brand choice among Pepsi, 7-Up and Coke | [`Q2_cola_mnl.m`](code/Q2_cola_mnl.m) |
| Q3 | Poisson regression for foot-and-mouth disease in Turkey | [`Q3_fmd_poisson.m`](code/Q3_fmd_poisson.m) |
| Q4 | Zero-inflated Poisson: excess zeros and likelihood | analytical (see report) |

Run everything with `run('HW2/code/run_all.m')`.

---

## Q1: Multinomial logit

With $U_{ij}=x_{ij}'\beta+\varepsilon_{ij}$ and i.i.d. Type-I extreme value errors, integrating over $\varepsilon_{ij}$ (substituting $t=e^{-\varepsilon_{ij}}$) gives

$$P_{ij}=\frac{\exp(x_{ij}'\beta)}{\sum_{k=1}^{J}\exp(x_{ik}'\beta)},\qquad L(\beta)=\prod_{i=1}^{N}\prod_{j=1}^{J}P_{ij}^{\,y_{ij}}.$$

For the three-alternative model with one covariate:

$$\frac{\partial P_{i1}}{\partial x_{i1,2}}=\beta_2P_{i1}(1-P_{i1}),\qquad \frac{\partial P_{i2}}{\partial x_{i1,2}}=-\beta_2P_{i1}P_{i2},\qquad \frac{\partial P_{i3}}{\partial x_{i1,2}}=-\beta_2P_{i1}P_{i3}.$$

## Q2: Soft-drink choice

Coke Classic is the base alternative. Prices enter relative to Coke's price.

| | Pepsi vs. Coke | 7-Up vs. Coke |
|---|---:|---:|
| Intercept | 0.2819 (t = 4.39) | 0.0898 (t = 1.29) |
| Price diff. (Pepsi − Coke) | −2.6100 (t = −10.27) | 0.1716 (t = 0.73) |
| Price diff. (7-Up − Coke) | 0.5512 (t = 2.47) | −2.0653 (t = −9.59) |

**Predicted shares** at Pepsi = $1.00, 7-Up = $1.25, Coke = $1.10: Pepsi **0.511**, Coke 0.273, 7-Up 0.216.

**Effect of a $0.10 price increase**

| Price raised | ΔP(Coke) | ΔP(Pepsi) | ΔP(7-Up) |
|---|---:|---:|---:|
| Pepsi | +0.0294 | **−0.0669** | +0.0351 |
| 7-Up | +0.0087 | +0.0319 | **−0.0382** |
| Coke | **−0.0381** | +0.0351 | +0.0031 |

Each brand loses share when its own price rises, and the lost share flows to its competitors.

## Q3: Foot-and-mouth disease, 1998

| Region | Variable | Mean | Median | SD | Min | Max |
|---|---|---:|---:|---:|---:|---:|
| Western | FMD cases | 1.075 | 0 | 1.654 | 0 | 5 |
| Western | Cattle | 157,850 | 135,025 | 82,619 | 39,330 | 392,840 |
| Eastern | FMD cases | 1.19 | 0 | 2.40 | 0 | 9 |
| Eastern | Cattle | 170,080 | 121,675 | 148,099 | 43,670 | 622,650 |

$$\log\lambda_i=\beta_1+\beta_2\ln(\text{Cattle}_i)+\beta_3\,\text{EasternTurkey}_i$$

| Variable | Coef. | t | 95% CI |
|---|---:|---:|---|
| Intercept | −11.254 | −4.94 | (−15.723, −6.785) |
| ln(Cattle) | 0.947 | 5.03 | (0.578, 1.316) |
| EasternTurkey | 0.041 | 0.17 | (−0.425, 0.510) |

Herd size matters: since the regressor is in logs, β₂ is an elasticity, so 1% more cattle means about 0.95% more expected outbreaks. Region has no significant effect once herd size is controlled for.

| Scenario | E[FMD] | Pr(FMD ≤ 2) |
|---|---:|---:|
| Eastern, Cattle = 170,080 | 1.22 | 0.875 |
| Western, Cattle = 157,850 | 1.09 | 0.902 |

## Q4: Zero-inflated Poisson

$$\Pr(y=0)=\pi+(1-\pi)e^{-\lambda},\qquad \Pr(y=0)_{\text{ZIP}}-\Pr(y=0)_{\text{Poisson}}=\pi(1-e^{-\lambda})>0$$

$$L(\pi,\lambda\mid y)=\prod_{y_i=0}\big[\pi+(1-\pi)e^{-\lambda}\big]\ \prod_{y_i>0}(1-\pi)\frac{\lambda^{y_i}e^{-\lambda}}{y_i!}$$
