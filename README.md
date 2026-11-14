<div align="center">

# ECO351 · Econometrics II

**Microeconometrics in MATLAB: binary & ordered choice, multinomial logit, count data, censoring & truncation, sample selection, treatment effects, and panel data.**

![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-0076A8?logo=mathworks&logoColor=white)
![Toolboxes](https://img.shields.io/badge/Toolboxes-Statistics%20%7C%20Optimization%20%7C%20Symbolic-informational)
![Homeworks](https://img.shields.io/badge/homeworks-4-success)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

*IIT Kanpur · Semester I, 2025–26 · Instructor: Prof. M. A. Rahman*

</div>

---

## Overview

This repository holds my solutions to the four homework assignments for **ECO351 (Econometrics II)**. Each homework mixes analytical derivations (likelihoods, marginal effects, estimators) with empirical work on real datasets. All estimation is done in MATLAB, as the course requires, and most of the likelihoods are coded by hand (Newton–Raphson, `fminunc`) instead of relying on canned routines.

| # | Topics | Datasets | Highlights |
|:-:|--------|----------|------------|
| [**HW1**](HW1/) | OLS vs. ML, binary probit & logit, ordinal probit | IgG, Pew marijuana survey, NLSY educational attainment | Symbolic Newton–Raphson MLE; average covariate effects; ordinal probit with reparameterised cut-points |
| [**HW2**](HW2/) | Multinomial / conditional logit, Poisson regression, zero-inflated Poisson | Cola brand choice, FMD outbreaks in Turkey | Gumbel → MNL derivation; cross-price effects; ZIP likelihood |
| [**HW3**](HW3/) | Truncation & censoring, Tobit, Heckman selection, endogenous treatment effects | Aptitude scores, Mroz (1987), NSW job training (jtrain2) | Hand-coded Tobit, Heckman two-step + bootstrap + FIML, bivariate-normal treatment model |
| [**HW4**](HW4/) | Selection-model marginal effects, fixed effects / within estimator, first differences | US state murder rates & executions | Within estimator derivation; HC1 and cluster-robust SEs; outlier (Texas) sensitivity |

## Repository layout

```
ECO351-Econometrics-II/
├── README.md               ← you are here
├── data/                   ← put the course .xlsx datasets here (not committed, see data/README.md)
├── HW1/
│   ├── README.md           ← questions, approach and key results
│   ├── assignment/         ← original problem set (PDF)
│   ├── report/             ← typeset solutions (PDF, LaTeX)
│   ├── code/               ← MATLAB scripts, one per question, plus run_all.m
│   └── figures/            ← plots produced by the code
├── HW2/ …
├── HW3/ …
└── HW4/ …
```

## Getting started

1. **Clone** the repository.
   ```bash
   git clone https://github.com/<your-username>/ECO351-Econometrics-II.git
   ```
2. **Add the datasets.** Copy the course data files into `data/`. The expected file names are listed in [`data/README.md`](data/README.md).
3. **Run a homework** from MATLAB:
   ```matlab
   run('HW3/code/run_all.m')      % or run any single script, e.g. HW3/code/Q3_heckman.m
   ```
   Every script finds `data/` from its own location, so you can run it from any working directory.

### Requirements

- MATLAB R2023b or newer (scripts use local functions at the end of script files)
- Statistics and Machine Learning Toolbox: `glmfit`, `fitglm`, `fitlm`, `mnrfit`, `normcdf`, …
- Optimization Toolbox: `fminunc`
- Symbolic Math Toolbox: only for HW1 Q1(d)

## Methods at a glance

<details>
<summary><b>Discrete choice</b></summary>

- Binary probit and logit by MLE, with average covariate effects for continuous and discrete changes (HW1)
- Ordinal probit with monotone cut-points `κ₁ = a₁, κⱼ = a₁ + Σ exp(aₖ)`, plus hit-rate and McFadden's R² (HW1)
- Multinomial logit derived from Type-I extreme value errors, with own- and cross-marginal effects `β Pᵢⱼ(1 − Pᵢⱼ)` and `−β Pᵢⱼ Pᵢₖ` (HW2)
</details>

<details>
<summary><b>Count data</b></summary>

- Poisson regression with log link, 95% Wald CIs, predicted means and `Pr(Y ≤ 2)` (HW2)
- Zero-inflated Poisson: excess-zero probability `π + (1 − π)e^{−λ}` and sample likelihood (HW2)
</details>

<details>
<summary><b>Limited dependent variables & selection</b></summary>

- Exponential duration model under truncation vs. censoring (HW3)
- Tobit Type I with right-censoring at 800, estimated with `fminunc` on `[β, log σ]` (HW3)
- Heckman sample selection: two-step with IMR, FGLS / robust SEs, a 400-draw full-sample bootstrap, and FIML with `ρ = tanh(·)` and `σ = exp(·)` (HW3)
- Endogenous binary treatment-effect model by FIML (HW3)
- Marginal effects in the selection model for hurdle-only, outcome-only and shared covariates (HW4)
</details>

<details>
<summary><b>Panel data</b></summary>

- Fixed-effects within estimator derived from scratch with the `I_n ⊗ M_T` demeaning operator (HW4)
- Pooled OLS vs. first differences vs. within estimation, with HC1 and state-clustered SEs (HW4)
</details>

## Authors

- **Himanshu Mahale** (230476)
- **Vineet Nagrale** (231158): group partner on HW2–HW4

## License

The code is released under the [MIT License](LICENSE). The problem sets in `*/assignment/` belong to the course instructor and are included for reference only. The datasets are not redistributed here.
