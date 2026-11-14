%% Q1_igg.m
% ECO351 - Homework 1, Question 1
% Serum immunoglobulin-G (IgG) vs. age for 298 children.
%   (a) descriptive statistics + histograms
%   (b) OLS:  IgG = b1 + b2*age + u
%   (c) OLS:  IgG = b1 + b2*age + b3*age^2 + e
%   (d) ML estimation of Model 2 with a hand-coded Newton-Raphson routine
%
% Data: data/igg.xlsx  (see data/README.md)

clear; clc; close all;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
outDir  = fullfile(here, '..', 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

T  = readtable(fullfile(dataDir, 'igg.xlsx'));
vn = lower(string(T.Properties.VariableNames));
idx_igg = find(vn == "igg", 1);
idx_age = find(vn == "age", 1);
if isempty(idx_igg) || isempty(idx_age)
    error('Could not find variables "igg" and/or "age" in igg.xlsx.');
end

IgG = T{:, idx_igg};
age = T{:, idx_age};
valid = ~isnan(IgG) & ~isnan(age);
IgG = IgG(valid); age = age(valid); n = numel(IgG);

disp('================================================================');
disp('===== Question 1: IgG vs Age (subparts a-d)');
disp('================================================================');

%% (a) Descriptive statistics and histograms
desc_stats = @(x) [mean(x), median(x), std(x), max(x), min(x)];
S = array2table([desc_stats(IgG); desc_stats(age)], ...
    'VariableNames', {'mean','median','sd','max','min'}, ...
    'RowNames',      {'IgG','age'});
disp(S);

skew_igg = skewness(IgG); skew_age = skewness(age);
fprintf('Skewness: IgG = %.4f ; age = %.4f\n', skew_igg, skew_age);

figure('Visible','off'); histogram(IgG, 20); title('Histogram: IgG'); xlabel('IgG'); ylabel('Count');
saveas(gcf, fullfile(outDir, 'q1a_hist_igg.png'));
figure('Visible','off'); histogram(age, 20); title('Histogram: Age'); xlabel('Age'); ylabel('Count');
saveas(gcf, fullfile(outDir, 'q1a_hist_age.png'));

%% (b) Model 1: linear regression
X1 = [ones(n,1), age];
[b1_vec, ~, r1, ~, stats1] = regress(IgG, X1);
sigma2_r1 = var(r1);
se1 = sqrt(diag(sigma2_r1 * inv(X1' * X1)));
t1  = b1_vec ./ se1;  R2_1 = stats1(1);
tbl1 = table({'Intercept';'Age'}, b1_vec, se1, t1, ...
    'VariableNames', {'term','coef','se','t'});
tbl1.R2 = [R2_1; R2_1];
disp(tbl1);

%% (c) Model 2: quadratic regression
X2 = [ones(n,1), age, age.^2];
[b2_vec, ~, r2, ~, stats2] = regress(IgG, X2);
sigma2_r2 = var(r2);
se2 = sqrt(diag(sigma2_r2 * inv(X2' * X2)));
t2  = b2_vec ./ se2;  R2_2 = stats2(1);
tbl2 = table({'Intercept';'Age';'Age^2'}, b2_vec, se2, t2, ...
    'VariableNames', {'term','coef','se','t'});
tbl2.R2 = [R2_2; R2_2; R2_2];
disp(tbl2);

%% (d) Maximum likelihood (errors ~ N(0,1)) via Newton-Raphson
syms b0 b1 b2 real
y_num = IgG(:); x_num = age(:); n_num = double(n);
mu_sym  = b0 + b1 .* x_num + b2 .* (x_num.^2);
nll_sym = 0.5 * sum((y_num - mu_sym).^2) + 0.5 * n_num * log(2*pi);   % negative log-likelihood

[grad_sym, hess_sym] = symbolic_gradient_hessian(nll_sym, [b0; b1; b2]);
Xinit = double(b2_vec);                                               % start from OLS
[xopt, ~, ~, ~] = symbolic_newton_raphson([b0; b1; b2], Xinit, grad_sym, hess_sym);

b_hat = double(xopt(:));
H_val = double(subs(hess_sym, [b0; b1; b2], xopt));
if rcond(H_val) < 1e-12, cov_b = pinv(H_val); else, cov_b = inv(H_val); end
se_b = sqrt(abs(diag(cov_b)));
t_b  = b_hat ./ se_b;

mle_tbl = table({'Intercept';'Age';'Age^2'}, b_hat, se_b, t_b, ...
    'VariableNames', {'term','mle','se','t'});
disp(mle_tbl);

%% ---------------------------------------------------------------------
%% Local functions
%% ---------------------------------------------------------------------
function [grad, hess] = symbolic_gradient_hessian(f, X)
    N = size(X,1);
    grad = sym('temp', [N 1]);
    for ii = 1:N, grad(ii) = simplify(diff(f, X(ii))); end
    hess = sym('temp', [N N]);
    for ii = 1:N
        for jj = 1:N
            hess(ii,jj) = simplify(diff(grad(ii), X(jj)));
        end
    end
end

function [xopt, num_iter, status_code, status_message] = symbolic_newton_raphson(X, Xinit, grad, hess)
    epsilon = 1e-15; max_iterations = 100; N = size(grad,1);
    xval = Xinit; err = inf; count_iterations = 0;
    gradval = double(subs(grad, X, xval));
    hessval = double(subs(hess, X, xval));

    while (err > epsilon) && (count_iterations <= max_iterations)
        if rank(hessval) == N
            step = hessval \ gradval;
        else
            hess_p = hessval(:, 1:rank(hessval));
            step = zeros(N,1);
            step(1:rank(hessval)) = hess_p \ gradval;
        end
        if any(isnan(step) | isinf(step)), break; end

        xval    = xval - step;
        gradval = double(subs(grad, X, xval));
        hessval = double(subs(hess, X, xval));
        err = (gradval' * gradval)^0.5;
        count_iterations = count_iterations + 1;
    end

    xopt = xval; num_iter = count_iterations;
    if err <= epsilon
        status_message = 'success';             status_code = 0;
    elseif count_iterations > max_iterations
        status_message = 'failure - max iters'; status_code = 1;
    else
        status_message = 'failure - hessian';   status_code = 2;
    end
end
