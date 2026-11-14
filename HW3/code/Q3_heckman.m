%% Q3_heckman.m
% ECO351 - Homework 3, Question 3
% Married women's hours of work (Mroz data) and sample-selection correction.
%   (a) OLS on all observations           (b) OLS on hours > 0 only
%   (c) probit for labour-force participation
%   (d) inverse Mills ratio from the probit
%   (e) Heckman two-step: OLS with IMR     (f) heteroscedasticity-robust (FGLS) SEs
%   (g) full-sample bootstrap SEs, B = 400 (h) Heckman selection model by FIML
%
% Data: data/mroz.xlsx  (sheet "Data")  -- see data/README.md

clear; clc; close all;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
mroz    = readtable(fullfile(dataDir, 'mroz.xlsx'), 'Sheet', 'Data');
n       = height(mroz);

%% (a) OLS on all observations
fprintf('\n===== (a) OLS, all observations =====\n');
model = fitlm(mroz, 'hours ~ educ + exper + hhours');
disp(model);
res_a = model.Residuals.Raw;
disp(table(["educ";"exper";"hhours"], ...
    [corr(res_a, mroz.educ); corr(res_a, mroz.exper); corr(res_a, mroz.hhours)], ...
    'VariableNames', {'Variable','Corr_with_Residuals'}));

%% (b) OLS on working women only (hours > 0)
fprintf('\n===== (b) OLS, hours > 0 =====\n');
Tpos   = mroz(mroz.hours > 0, :);
model1 = fitlm(Tpos, 'hours ~ educ + exper + hhours');
disp(model1);
res_b = model1.Residuals.Raw;
disp(table(["educ";"exper";"hhours"], ...
    [corr(res_b, Tpos.educ); corr(res_b, Tpos.exper); corr(res_b, Tpos.hhours)], ...
    'VariableNames', {'Variable','Corr_with_Residuals'}));

%% (c) Probit for labour-force participation
fprintf('\n===== (c) Probit: LFP =====\n');
lfp = mroz.lfp;
X1  = [mroz.exper, mroz.kidsl6, mroz.kids618, mroz.mtr, mroz.largecity];
[beta_probit, ~, stats_probit] = glmfit(X1, lfp, 'binomial', 'Link', 'probit');
disp(table({'Intercept';'exper';'kidsl6';'kids618';'MTR';'largecity'}, ...
    beta_probit, stats_probit.se, stats_probit.t, ...
    'VariableNames', {'Variable','Coef','SE','t'}));

%% (d) Inverse Mills ratio
w      = [ones(n,1), X1] * beta_probit;
lambda = normpdf(w) ./ normcdf(w);
fprintf('\n===== (d) Inverse Mills ratio =====\n');
fprintf('mean(lambda) = %.4f,  var(lambda) = %.4f\n', mean(lambda), var(lambda));

%% (e) Heckman two-step: OLS with the IMR on hours > 0
fprintf('\n===== (e) Two-step: OLS with IMR =====\n');
idx        = mroz.hours > 0;
hours_pos  = mroz.hours(idx);
lambda_pos = lambda(idx);
T2 = table(hours_pos, mroz.educ(idx), mroz.exper(idx), mroz.hhours(idx), lambda_pos, ...
    'VariableNames', {'hours','educ','exper','hhours','lambda'});
model2 = fitlm(T2, 'hours ~ educ + exper + hhours + lambda');
disp(model2);

X2 = [ones(sum(idx),1), T2.educ, T2.exper, T2.hhours, T2.lambda];
beta2 = regress(hours_pos, X2);
e   = hours_pos - X2 * beta2;
[nn, k] = size(X2);
s2  = (e' * e) / (nn - k);
se_beta2 = sqrt(diag(s2 * inv(X2' * X2)));
disp(table(beta2, se_beta2, beta2 ./ se_beta2, 'VariableNames', {'Estimate','Std_Error','t_Stat'}));

%% (f) Heteroscedasticity-robust (FGLS) estimates
fprintf('\n===== (f) FGLS with robust standard errors =====\n');
y1    = hours_pos;
b_ols = (X2' * X2) \ (X2' * y1);
u_hat = y1 - X2 * b_ols;
Omega_hat  = diag(u_hat.^2);
b_fgls     = (X2' * (Omega_hat \ X2)) \ (X2' * (Omega_hat \ y1));
Var_b_fgls = inv(X2' * X2) * (X2' * Omega_hat * X2) * inv(X2' * X2);
se_fgls    = sqrt(diag(Var_b_fgls));
disp(table(b_fgls, se_fgls, b_fgls ./ se_fgls, 'VariableNames', {'Estimate','Robust_SE','t_Stat'}));

%% (g) Full-sample bootstrap of the two-step estimator (B = 400)
fprintf('\n===== (g) Bootstrap standard errors (B = 400) =====\n');
X1_full  = [mroz.exper, mroz.kidsl6, mroz.kids618, mroz.mtr, mroz.largecity];
lfp_full = double(mroz.hours > 0);
lambda_full = imr(X1_full, glmfit(X1_full, lfp_full, 'binomial', 'link', 'probit'));

idx_work = mroz.hours > 0;
y_w  = mroz.hours(idx_work);
X2_w = [ones(sum(idx_work),1), mroz.educ(idx_work), mroz.exper(idx_work), ...
        mroz.hhours(idx_work), lambda_full(idx_work)];
beta2_hat = (X2_w' * X2_w) \ (X2_w' * y_w);
k2 = size(X2_w, 2);
resid_w = y_w - X2_w * beta2_hat;
se_ols  = sqrt(diag((resid_w' * resid_w) / (sum(idx_work) - k2) * inv(X2_w' * X2_w)));

B = 400;
boot_b = nan(B, k2);
rng(123);                                       % reproducible
for b = 1:B
    draw   = randsample(n, n, true);            % resample women with replacement
    mroz_b = mroz(draw, :);
    X1_b   = [mroz_b.exper, mroz_b.kidsl6, mroz_b.kids618, mroz_b.mtr, mroz_b.largecity];
    lfp_b  = double(mroz_b.hours > 0);
    try
        g_b = glmfit(X1_b, lfp_b, 'binomial', 'link', 'probit');
    catch
        continue;                               % skip draw if probit fails
    end
    lambda_b = imr(X1_b, g_b);
    idx_b = mroz_b.hours > 0;
    if sum(idx_b) < 5, continue; end
    X2_b = [ones(sum(idx_b),1), mroz_b.educ(idx_b), mroz_b.exper(idx_b), ...
            mroz_b.hhours(idx_b), lambda_b(idx_b)];
    boot_b(b,:) = ((X2_b' * X2_b) \ (X2_b' * mroz_b.hours(idx_b)))';
end

boot_good = boot_b(~any(isnan(boot_b), 2), :);
B_eff = size(boot_good, 1);
if B_eff == 0, error('No successful bootstrap draws.'); end
if B_eff < 30, warning('Only %d successful bootstrap draws.', B_eff); end

boot_se   = std(boot_good, 0, 1)';
boot_bias = mean(boot_good, 1)' - beta2_hat;
t_boot    = beta2_hat ./ boot_se;
p_norm    = 2 * (1 - normcdf(abs(t_boot)));
ci        = prctile(boot_good, [2.5, 97.5]);

varnames = {'Intercept','educ','exper','hhours','lambda'};
disp(table(varnames', beta2_hat, se_ols, boot_se, t_boot, p_norm, ci(1,:)', ci(2,:)', boot_bias, ...
    'VariableNames', {'Variable','Estimate','OLS_SE','Bootstrap_SE','t_boot','p_norm', ...
                      'CI_lower_95','CI_upper_95','Boot_Bias'}));
fprintf('Requested B = %d, successful draws = %d\n', B, B_eff);

%% (h) Heckman selection model by FIML
fprintf('\n===== (h) Heckman selection model (FIML) =====\n');
y   = mroz.hours;
lfp = double(mroz.hours > 0);
X = [ones(n,1), mroz.educ, mroz.exper, mroz.hhours];                                 % outcome eq.
Z = [ones(n,1), mroz.exper, mroz.kidsl6, mroz.kids618, mroz.mtr, mroz.largecity];    % selection eq.
kX = size(X,2);  kZ = size(Z,2);

% Starting values from the two-step estimator
gamma_init  = glmfit(Z(:,2:end), lfp, 'binomial', 'link', 'probit');
lambda_init = imr(Z(:,2:end), gamma_init);
sel         = lfp == 1;
X2_work     = [X(sel,:), lambda_init(sel)];
b_2s        = (X2_work' * X2_work) \ (X2_work' * y(sel));
beta_init   = b_2s(1:kX);
resid_init  = y(sel) - X2_work * b_2s;
sigma_init  = sqrt(sum(resid_init.^2) / max(1, sum(sel) - size(X2_work,2)));
rho_init    = b_2s(end) / sigma_init;           % beta_lambda = rho * sigma
if abs(rho_init) >= 0.98, rho_init = sign(rho_init) * 0.5; end

theta0 = [beta_init; gamma_init; atanh(rho_init); log(sigma_init)];

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', ...
    'MaxFunctionEvaluations', 2e5, 'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8);
obj = @(th) heckman_nll(th, y, X, lfp, Z);
[theta_hat, neglogL, exitflag, output, ~, hessian] = fminunc(obj, theta0, options);
if exitflag <= 0
    warning('fminunc did not converge (exitflag = %d).', exitflag);
end

rcond_h = rcond(hessian);
if rcond_h < 1e-12
    warning('Hessian poorly conditioned (rcond = %.3e); using pseudo-inverse.', rcond_h);
    vcov_theta = pinv(hessian);
else
    vcov_theta = inv(hessian);
end

pos_beta  = 1:kX;
pos_gamma = kX + (1:kZ);
pos_rho   = kX + kZ + 1;
pos_sig   = numel(theta_hat);

beta_hat  = theta_hat(pos_beta);
gamma_hat = theta_hat(pos_gamma);
rho_hat   = tanh(theta_hat(pos_rho));
sigma_hat = exp(theta_hat(pos_sig));

se_beta  = sqrt(diag(vcov_theta(pos_beta,  pos_beta)));
se_gamma = sqrt(diag(vcov_theta(pos_gamma, pos_gamma)));
se_rho   = sqrt(max((1 - rho_hat^2)^2 * vcov_theta(pos_rho, pos_rho), 0));    % delta method
se_sigma = sqrt(max(sigma_hat^2      * vcov_theta(pos_sig, pos_sig), 0));

report = @(names, est, se) table(names', est, se, est./se, 2*(1-normcdf(abs(est./se))), ...
    'VariableNames', {'Variable','Estimate','StdError','tStat','pValue'});

fprintf('\n--- FIML: OUTCOME EQUATION ---\n');
disp(report({'Intercept','educ','exper','hhours'}, beta_hat, se_beta));
fprintf('\n--- FIML: SELECTION EQUATION ---\n');
disp(report({'Intercept','exper','kidsl6','kids618','MTR','largecity'}, gamma_hat, se_gamma));
fprintf('\n--- CORRELATION STRUCTURE ---\n');
disp(report({'rho','sigma'}, [rho_hat; sigma_hat], [se_rho; se_sigma]));

fprintf('\nLog-likelihood = %.6f\n', -neglogL);
fprintf('fminunc exitflag = %d, iterations = %d, funcCount = %d\n', ...
    exitflag, output.iterations, output.funcCount);
fprintf('Hessian rcond = %.3e\n', rcond_h);

%% ---------------------------------------------------------------------
%% Local functions
%% ---------------------------------------------------------------------
function lam = imr(Xnoint, g)
    % Inverse Mills ratio phi(w)/Phi(w) for probit index w = [1 X]*g
    w   = [ones(size(Xnoint,1),1), Xnoint] * g;
    Phi = max(min(normcdf(w), 1 - 1e-10), 1e-10);
    lam = normpdf(w) ./ Phi;
end

function nll = heckman_nll(theta, y, X, lfp, Z)
    % theta = [beta (kX); gamma (kZ); atanh(rho); log(sigma)]
    kX = size(X,2);  kZ = size(Z,2);
    beta  = theta(1:kX);
    gamma = theta(kX+1:kX+kZ);
    rho   = tanh(theta(kX+kZ+1));
    sigma = exp(theta(end));
    if sigma <= 0 || abs(rho) >= 0.9999999
        nll = 1e20; return;
    end

    xb = X * beta;
    zg = Z * gamma;
    sel = (lfp == 1);

    % Working women: (1/sigma) phi(eps) Phi((z'g + rho*eps)/sqrt(1-rho^2))
    eps_obs = (y(sel) - xb(sel)) / sigma;
    v_obs   = (zg(sel) + rho .* eps_obs) ./ sqrt(1 - rho^2);
    logL_obs = log((1/sigma) .* normpdf(eps_obs) .* max(normcdf(v_obs), 1e-12));

    % Non-working women: Pr(LFP = 0) = Phi(-z'g)
    logL_unobs = log(max(normcdf(-zg(~sel)), 1e-12));

    nll = -(sum(logL_obs) + sum(logL_unobs));
end
