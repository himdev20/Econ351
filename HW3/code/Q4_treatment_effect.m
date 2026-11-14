%% Q4_treatment_effect.m
% ECO351 - Homework 3, Question 4
% Effect of job training on 1978 earnings (jtrain2, n = 445).
%   (a) fraction trained            (b) mean re78 by training status
%   (c) unemployment rate by status (d) comment on causality
%   (f) endogenous treatment-effect model estimated by FIML
%         T* = Z'gamma + u1,  T = 1(T* > 0)
%         re78 = X'beta + alpha*T + u2,  (u1,u2) ~ BVN(0, [1 rho*s; rho*s s^2])
%   (g) inference on alpha
%
% Data: data/jtrain2.xlsx  (sheet "Data")  -- see data/README.md

clear; clc; close all;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
T = readtable(fullfile(dataDir, 'jtrain2.xlsx'), 'Sheet', 'Data');

train    = T.train;     re78 = T.re78;   re74 = T.re74;   re75 = T.re75;
age      = T.age;       educ = T.educ;   black = T.black; hisp = T.hisp;
married  = T.married;   nodegree = T.nodegree;
unemp74  = T.unem74;    unem78 = T.unem78;
n = height(T);

%% (a)-(c) Descriptive comparisons
fprintf('(a) Fraction receiving training = %.4f  (%d / %d)\n', mean(train), sum(train), n);
fprintf('(b) Mean re78: trained = %.4f ; untrained = %.4f\n', ...
    mean(re78(train==1)), mean(re78(train==0)));
fprintf('(c) Fraction unemployed in 1978: trained = %.4f ; untrained = %.4f\n', ...
    mean(unem78(train==1)), mean(unem78(train==0)));

%% (d) Comment
fprintf(['(d) The raw differences in (b) and (c) favour the trained group, but they are not\n' ...
         '    causal: participation was not randomised, so selection on unobservables\n' ...
         '    must be modelled (treatment-effect model below) or ruled out.\n\n']);

%% (f) Treatment-effect model by FIML
Z = [ones(n,1), age, educ, black, hisp, married, nodegree, re74, unemp74];   % treatment eq.
X = [ones(n,1), age, educ, black, hisp, married, re74, re75];                % outcome eq.
kZ = size(Z,2);  kX = size(X,2);

% Starting values: probit for gamma, OLS for (beta, alpha), rho = 0
gamma0 = glmfit(Z(:,2:end), train, 'binomial', 'link', 'probit');
b_tmp  = [X, train] \ re78;
beta0  = b_tmp(1:kX);  alpha0 = b_tmp(end);
sigma0 = std(re78 - X*beta0 - alpha0*train);
theta0 = [gamma0; beta0; alpha0; atanh(0); log(sigma0)];

negloglik = @(theta) negloglik_treatment(theta, re78, X, train, Z);
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', ...
    'MaxFunctionEvaluations', 5e4, 'MaxIterations', 1e4);
[theta_hat, fval, exitflag, output, grad, hess] = fminunc(negloglik, theta0, options);

gamma_hat = theta_hat(1:kZ);
beta_hat  = theta_hat(kZ+1:kZ+kX);
alpha_hat = theta_hat(kZ+kX+1);
idx_rho   = kZ + kX + 2;
idx_sig   = numel(theta_hat);
rho_hat   = tanh(theta_hat(idx_rho));
sigma_hat = exp(theta_hat(idx_sig));

vcov_theta = inv(hess);
se_theta   = sqrt(diag(vcov_theta));
se_rho     = abs(se_theta(idx_rho)) * sech(theta_hat(idx_rho))^2;   % delta method
se_sigma   = abs(se_theta(idx_sig)) * sigma_hat;

fprintf('\n--- FIML Treatment Model Estimates ---\n');
zNames = {'Intercept','age','educ','black','hisp','married','nodegree','re74','unemp74'};
xNames = {'Intercept','age','educ','black','hisp','married','re74','re75'};
se_g = se_theta(1:kZ);  se_b = se_theta(kZ+1:kZ+kX);  se_a = se_theta(kZ+kX+1);

fprintf('\nTreatment (selection) equation:\n');
disp(table(zNames', gamma_hat, se_g, gamma_hat./se_g, 'VariableNames', {'Variable','Estimate','StdError','tStat'}));
fprintf('Outcome equation (re78):\n');
disp(table([xNames, {'alpha (train)'}]', [beta_hat; alpha_hat], [se_b; se_a], [beta_hat; alpha_hat]./[se_b; se_a], ...
    'VariableNames', {'Variable','Estimate','StdError','tStat'}));

fprintf('Covariance components:\n');
fprintf('  rho        = %8.4f  (se = %.4f, t = %.3f)\n', rho_hat, se_rho, rho_hat/se_rho);
fprintf('  sigma      = %8.4f  (se = %.4f, t = %.3f)\n', sigma_hat, se_sigma, sigma_hat/se_sigma);
fprintf('  sigma_12   = %8.4f  (= rho * sigma)\n', rho_hat * sigma_hat);

%% (g) Treatment effect
t_alpha = alpha_hat / se_a;
fprintf('\n(g) alpha = %.4f, se = %.4f, t = %.3f\n', alpha_hat, se_a, t_alpha);
if abs(t_alpha) > 1.96
    fprintf('    => alpha is statistically significant at the 5%% level\n');
else
    fprintf('    => alpha is NOT statistically significant at the 5%% level\n');
end

%% ---------------------------------------------------------------------
%% Local functions
%% ---------------------------------------------------------------------
function negL = negloglik_treatment(theta, y, X, T, Z)
    % theta = [gamma (kZ); beta (kX); alpha; atanh(rho); log(sigma)]
    kZ = size(Z,2);  kX = size(X,2);
    gamma = theta(1:kZ);
    beta  = theta(kZ+1:kZ+kX);
    alpha = theta(kZ+kX+1);
    rho   = tanh(theta(kZ+kX+2));
    sigma = exp(theta(end));

    z  = Z * gamma;
    xb = X * beta;
    sel = (T == 1);
    tiny = 1e-200;

    % Treated: f(y | T=1) * Pr(T=1 | u2)
    eps1 = (y(sel) - xb(sel) - alpha) ./ sigma;
    c1   = normpdf(eps1) ./ sigma .* normcdf((z(sel) + rho .* eps1) ./ sqrt(1 - rho^2));
    c1(c1 < tiny) = tiny;

    % Untreated: f(y | T=0) * Pr(T=0 | u2)
    eps0 = (y(~sel) - xb(~sel)) ./ sigma;
    c0   = normpdf(eps0) ./ sigma .* normcdf((-z(~sel) - rho .* eps0) ./ sqrt(1 - rho^2));
    c0(c0 < tiny) = tiny;

    negL = -(sum(log(c1)) + sum(log(c0)));
end
