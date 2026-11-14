%% Q2_tobit.m
% ECO351 - Homework 3, Question 2
% Tobit Type I (censored regression) for academic aptitude, right-censored at 800.
%   (a) program-type shares            (b) descriptive summary
%   (c) histograms of apt by program   (d) correlations + scatter plots
%   (e) Tobit MLE (academic = base)    (f) average partial effects on E[apt|x]
%   (g) residual and actual-vs-fitted plots
%
% Data: data/Aptitude.xlsx  (see data/README.md)

clear; clc; close all;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
figDir  = fullfile(here, '..', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

T = readtable(fullfile(dataDir, 'Aptitude.xlsx'));
T.Properties.VariableNames = lower(T.Properties.VariableNames);

prog  = string(T.prog);
n_row = height(T);
read_num = T.read;
math_num = T.math;
apt_num  = T.apt;                         % dependent variable

%% (a) Program categories and shares
% The prompt codes general = 2; 0/1 dummies are used here because they are
% what the Tobit model in part (e) needs.
T.academic   = double(prog == "academic");
T.general    = double(prog == "general");
T.vocational = double(prog == "vocational");

counts = [sum(T.academic), sum(T.general), sum(T.vocational)];
shares = counts / n_row * 100;
disp('==================== Q2 (a): Program Shares ====================');
labels = {'Academic','General','Vocational'};
for k = 1:3
    fprintf('%-12s %5.0f   %6.2f %%\n', labels{k}, counts(k), shares(k));
end

%% (b) Descriptive summary
disp(' '); disp('==================== Q2 (b): Descriptive Summary ====================');
vars = {'read','math','apt'};
for k = 1:numel(vars)
    v = T.(vars{k});
    fprintf('%-5s mean=%.2f, median=%.2f, std=%.2f, max=%.0f, min=%.0f\n', ...
        vars{k}, mean(v), median(v), std(v), max(v), min(v));
end

%% (c) Histograms of apt by program type
progs = ["academic","general","vocational"];
for k = 1:3
    figure('Visible','off');
    histogram(apt_num(prog == progs(k)));
    title("Aptitude Scores - " + labels{k}); xlabel('Aptitude Score'); ylabel('Frequency');
    saveas(gcf, fullfile(figDir, "q2c_hist_apt_" + progs(k) + ".png"));
end

%% (d) Correlation matrix and scatter plots
disp(' '); disp('==================== Q2 (d): Correlation Matrix ====================');
CorrelationMatrix = corr([read_num, math_num, apt_num], 'rows', 'complete');
disp(array2table(CorrelationMatrix, 'VariableNames', vars, 'RowNames', vars));

pairs = {read_num, math_num, 'Read Score', 'Math Score', 'Read vs Math', 'q2d_scatter_read_math';
         read_num, apt_num,  'Read Score', 'Apt Score',  'Read vs Apt',  'q2d_scatter_read_apt';
         math_num, apt_num,  'Math Score', 'Apt Score',  'Math vs Apt',  'q2d_scatter_math_apt'};
for k = 1:size(pairs,1)
    figure('Visible','off');
    scatter(pairs{k,1}, pairs{k,2}, 'filled');
    xlabel(pairs{k,3}); ylabel(pairs{k,4}); title(pairs{k,5});
    saveas(gcf, fullfile(figDir, [pairs{k,6} '.png']));
end
% Censoring shows up as a horizontal band of points at apt = 800.

%% (e) Tobit Type I estimation
disp(' '); disp('==================== Q2 (e): Tobit Type I Model ====================');
X = [ones(n_row,1), read_num, math_num, T.general, T.vocational];   % academic = base
Censor_Point = 800;
d = (apt_num < Censor_Point);                   % 1 = uncensored, 0 = censored at 800

num_predictors = size(X,2);
b0 = X(d,:) \ apt_num(d);                       % OLS on uncensored obs as starting values
s0 = std(apt_num(d) - X(d,:) * b0);
p0 = [b0; log(max(s0, 1e-6))];                  % [beta; log(sigma)]

opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'final');
[param_estimates, NegLogL, exitflag, output, grad, Hessian] = fminunc( ...
    @(p) tobit_loglikelihood(p, apt_num, X, d, Censor_Point), p0, opts);

Beta_Estimates  = param_estimates(1:num_predictors);
Sigma_Estimate  = exp(param_estimates(end));
Standard_Errors = sqrt(diag(inv(Hessian)));
T_Statistics    = Beta_Estimates ./ Standard_Errors(1:num_predictors);
P_Values        = 2 * (1 - normcdf(abs(T_Statistics)));

Coefficient_Names = {'Intercept'; 'read'; 'math'; 'prog_general'; 'prog_vocational'};
Results = table(Beta_Estimates, Standard_Errors(1:num_predictors), T_Statistics, P_Values, ...
    'VariableNames', {'Coefficient','StdError','tStatistic','pValue'}, ...
    'RowNames', Coefficient_Names);
fprintf('Sigma: %.4f\n', Sigma_Estimate);
disp(Results);
fprintf('Are all the coefficients significant (p < 0.05)? %s\n', string(all(P_Values < 0.05)));

%% (f) Average partial effects on the observed outcome
disp(' '); disp('==================== Q2 (f): Average Partial Effects ====================');
% For a continuous regressor under right-censoring at C:
%   dE[y|x]/dx_k = b_k * Phi((C - x'b)/s)
mu = X * Beta_Estimates;
Phi_term = normcdf((Censor_Point - mu) / Sigma_Estimate);
APE_read = mean(Beta_Estimates(2) * Phi_term);
APE_math = mean(Beta_Estimates(3) * Phi_term);
fprintf('Expected change in apt for +1 READ: %.4f\n', APE_read);
fprintf('Expected change in apt for +1 MATH: %.4f\n', APE_math);

E_apt = @(Xm) expected_apt(Xm, Beta_Estimates, Sigma_Estimate, Censor_Point);
X_base = X; X_base(:,4:5) = 0;                  % everyone academic
X_voc  = X; X_voc(:,4) = 0; X_voc(:,5) = 1;     % everyone vocational
APE_vocational = mean(E_apt(X_voc) - E_apt(X_base));
fprintf('Expected change in apt for VOCATIONAL (vs. academic): %.4f\n', APE_vocational);

%% (g) Residuals vs fitted, actual vs fitted
pred_observed = E_apt(X);
residuals     = apt_num - pred_observed;

figure('Visible','off');
scatter(pred_observed, residuals, 'filled');
xlabel('Fitted apt  E[apt|X]'); ylabel('Residual  (apt - E[apt|X])');
title('Residuals vs Fitted Apt Values');
saveas(gcf, fullfile(figDir, 'q2g_residuals_vs_fitted.png'));

figure('Visible','off');
scatter(pred_observed, apt_num, 'filled');
xlabel('Fitted apt  E[apt|X]'); ylabel('Actual apt');
title('Actual vs Fitted Apt Values');
saveas(gcf, fullfile(figDir, 'q2g_actual_vs_fitted.png'));

%% ---------------------------------------------------------------------
%% Local functions
%% ---------------------------------------------------------------------
function neglogL = tobit_loglikelihood(p, Y, X, d, C)
    % Right-censored Tobit: p = [beta; log(sigma)]
    beta  = p(1:end-1);
    sigma = exp(p(end));
    mu    = X * beta;
    L_uncensored = log(normpdf(Y(d), mu(d), sigma));          % density for y < C
    L_censored   = log(1 - normcdf((C - mu(~d)) / sigma));    % Pr(y* >= C)
    neglogL = -sum(L_uncensored) - sum(L_censored);
end

function Ey = expected_apt(Xm, beta, sigma, C)
    % Expected observed aptitude used for the fitted values in parts (f)-(g)
    xb = Xm * beta;
    Ey = xb .* normcdf(xb / sigma) + sigma .* normpdf(xb / sigma) + C .* (1 - normcdf(xb / sigma));
end
