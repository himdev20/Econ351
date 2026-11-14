%% Q5_ordinal_probit.m
% ECO351 - Homework 1, Question 5
% Ordinal probit for educational attainment (4 ordered levels).
%   (a) Pr(y_i = j) = Phi(kappa_j - x_i'b) - Phi(kappa_{j-1} - x_i'b)
%       L(b, kappa) = prod_i prod_j Pr(y_i = j)^{1(y_i = j)}
%   (b) descriptive summary (mean/sd for continuous, count/% for categorical)
%   (c) MLE via numerical Newton-Raphson, hit-rate and McFadden's R^2
%   (d) average covariate effect of a $10,000 increase in family income
%
% Cut-points are re-parameterised as kappa_1 = a_1, kappa_j = a_1 + sum exp(a_2..a_j)
% so that they stay strictly increasing during optimisation.
%
% Data: data/EducationalAttainment.xlsx  (see data/README.md)

clear; clc; close all;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');

T = readtable(fullfile(dataDir, 'EducationalAttainment.xlsx'));
req = {'educLevel','famIncome','educMother','educFather','motherWorked','female', ...
       'black','urban','south','age15','age16','age17'};
for i = 1:numel(req)
    if ~ismember(req{i}, T.Properties.VariableNames)
        error('Missing variable %s in sheet.', req{i});
    end
end

y = T.educLevel;
famIncome = T.famIncome;
X = [sqrt(famIncome), T.educMother, T.educFather, T.motherWorked, ...
     T.female, T.black, T.urban, T.south, T.age15, T.age16, T.age17];
n = size(X,1);  K = numel(unique(y));

%% (a) Probability expression and likelihood -- see header comments
disp('--- (a) Ordinal probit probability and likelihood: see header ---');

%% (b) Descriptive summary
disp('--- (b) Descriptive summary ---');
cont_names = {'sqrt_famIncome','educMother','educFather'};
cont_data  = [sqrt(T.famIncome), T.educMother, T.educFather];
means = mean(cont_data, 'omitnan');  sds = std(cont_data, 'omitnan');

cat_names = {'motherWorked','female','black','urban','south','age15','age16','age17'};
cat_data  = [T.motherWorked, T.female, T.black, T.urban, T.south, T.age15, T.age16, T.age17];

rows = {};
for i = 1:numel(cont_names)
    rows(end+1,:) = {cont_names{i}, sprintf('mean=%.3f', means(i)), sprintf('sd=%.3f', sds(i))}; %#ok<SAGROW>
end
for i = 1:numel(cat_names)
    cnt = sum(cat_data(:,i) == 1);  pct = 100 * cnt / n;
    rows(end+1,:) = {cat_names{i}, sprintf('count=%d', cnt), sprintf('pct=%.2f%%', pct)}; %#ok<SAGROW>
end
edu_counts = accumarray(y, 1, [4,1]);  edu_pcts = 100 * edu_counts / n;
rows(end+1,:) = {'educLevel (1..4)', sprintf('counts=%s', mat2str(edu_counts')), ...
                 sprintf('pct=%s', mat2str(round(edu_pcts', 2)))};
DescTable = cell2table(rows, 'VariableNames', {'variable','mean_or_count','sd_or_pct'});
disp(DescTable);

%% (c) Ordinal probit by MLE
disp(' '); disp('--- (c) Estimate ordinal probit by MLE ---');
X_with_int = [ones(n,1) X];
p2 = size(X_with_int, 2);
theta0 = [zeros(p2,1); zeros(K-1,1)];

negloglik_fun = @(th) negloglik_ordprobit(y, X_with_int, th);
[theta_hat, H] = newton_raphson_numeric(negloglik_fun, theta0, 1e-8, 400);

beta_hat  = theta_hat(1:p2);
a_hat     = theta_hat(p2+1:end);
kappa_hat = a_to_kappa(a_hat);

cov_theta = inv(H);
se_theta  = sqrt(diag(cov_theta));
se_beta   = se_theta(1:p2);
t_beta    = beta_hat ./ se_beta;

names_out = {'Intercept','sqrt_famIncome','educMother','educFather','motherWorked', ...
             'female','black','urban','south','age15','age16','age17'};
coef_tbl = table(names_out', round(beta_hat,3), round(se_beta,3), round(t_beta,3), ...
    'VariableNames', {'term','coef','se','t'});
disp('Coefficient estimates:'); disp(coef_tbl);

[pr_matrix, yhat] = ordinal_probs_and_pred(X_with_int, beta_hat, kappa_hat);
hit_rate = mean(yhat == y);
LL_model = sum(log(max(pr_matrix(sub2ind(size(pr_matrix), (1:n)', y)), 1e-12)));
shares   = accumarray(y, 1) / n;
LL_null  = sum(log(max(shares(y), 1e-12)));
McFaddenR2 = 1 - LL_model / LL_null;
fprintf('Hit rate      = %.3f\n', hit_rate);
fprintf('McFadden R^2  = %.4f\n', McFaddenR2);

%% (d) Average covariate effect of +$10k family income (famIncome is in $1000s)
d  = 10;
X1 = X_with_int;  X2 = X_with_int;
X2(:,2) = sqrt(famIncome + d);
pr0 = ordinal_probs_and_pred(X1, beta_hat, kappa_hat);
pr1 = ordinal_probs_and_pred(X2, beta_hat, kappa_hat);
avg_dP = mean(pr1 - pr0, 1);
fprintf('\n(d) Average change in probabilities for +$10k famIncome:\n');
for j = 1:K
    fprintf('    P(y=%d): delta = %.4f\n', j, round(avg_dP(j), 4));
end

%% ---------------------------------------------------------------------
%% Local functions
%% ---------------------------------------------------------------------
function nll = negloglik_ordprobit(y, X, theta)
    p = size(X,2);
    beta = theta(1:p);  a = theta(p+1:end);
    kappa = a_to_kappa(a);  eta = X * beta;
    Kloc = numel(kappa) + 1;
    P = zeros(length(y), Kloc);
    P(:,1) = normcdf(kappa(1) - eta);
    for jj = 2:(Kloc-1)
        P(:,jj) = normcdf(kappa(jj) - eta) - normcdf(kappa(jj-1) - eta);
    end
    P(:,Kloc) = 1 - normcdf(kappa(end) - eta);
    idx  = sub2ind(size(P), (1:length(y))', y);
    pr_y = max(P(idx), 1e-12);
    nll  = -sum(log(pr_y));
end

function kappa = a_to_kappa(a)
    m = numel(a);  kappa = zeros(m,1);
    if m == 0, return; end
    kappa(1) = a(1);
    for j = 2:m, kappa(j) = a(1) + sum(exp(a(2:j))); end
end

function g = numgrad(f, theta)
    h = 1e-6;  g = zeros(numel(theta),1);  f0 = f(theta);
    for i = 1:numel(theta)
        th = theta;  th(i) = th(i) + h;
        g(i) = (f(th) - f0) / h;
    end
end

function H = numhess(f, theta)
    h = 1e-5;  k = numel(theta);  H = zeros(k);  f0 = f(theta);
    for i = 1:k
        ei = zeros(k,1); ei(i) = h;
        for j = i:k
            ej = zeros(k,1); ej(j) = h;
            H(i,j) = (f(theta+ei+ej) - f(theta+ei) - f(theta+ej) + f0) / h^2;
            H(j,i) = H(i,j);
        end
    end
end

function [theta_hat, H] = newton_raphson_numeric(fhandle, theta0, tol, maxit)
    theta = theta0(:);
    for it = 1:maxit
        g = numgrad(fhandle, theta);
        H = numhess(fhandle, theta);
        step = -H \ g;
        theta_new = theta + step;
        if norm(step) < tol * (1 + norm(theta)), theta = theta_new; break; end
        theta = theta_new;
    end
    theta_hat = theta;
    H = numhess(fhandle, theta_hat);
end

function [pr, yhat] = ordinal_probs_and_pred(X, beta, kappa)
    eta  = X * beta;  Kloc = numel(kappa) + 1;
    P = zeros(size(X,1), Kloc);
    P(:,1) = normcdf(kappa(1) - eta);
    for jj = 2:(Kloc-1)
        P(:,jj) = normcdf(kappa(jj) - eta) - normcdf(kappa(jj-1) - eta);
    end
    P(:,Kloc) = 1 - normcdf(kappa(end) - eta);
    pr = P;
    [~, yhat] = max(P, [], 2);
end
