%% Q3_fmd_poisson.m
% ECO351 - Homework 2, Question 3
% Poisson regression for foot-and-mouth disease (FMD) incidents in Turkish provinces, 1998.
%   log(lambda_i) = b1 + b2*ln(Cattle_i) + b3*EasternTurkey_i
%   (a) descriptive summary by region
%   (b) MLE, standard errors, t-values, 95% CIs
%   (d)-(g) expected counts and Pr(Y <= 2) for given cattle populations
%
% Columns of FMD.xlsx: [EasternTurkey, FMD1998, Cattle]
% Data: data/FMD.xlsx  (see data/README.md)

clc; clear;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
FMD     = readmatrix(fullfile(dataDir, 'FMD.xlsx'));

east    = FMD(:,1) == 1;
Eastern = FMD(east, :);
Western = FMD(~east, :);
fprintf('Provinces: %d Western, %d Eastern\n', size(Western,1), size(Eastern,1));

%% (a) Descriptive summary by region
summ = @(v) [mean(v), median(v), std(v), min(v), max(v)];
rowNames = {'Western FMD1998'; 'Western Cattle'; 'Eastern FMD1998'; 'Eastern Cattle'};
S = [summ(Western(:,2)); summ(Western(:,3)); summ(Eastern(:,2)); summ(Eastern(:,3))];
disp(array2table(S, 'VariableNames', {'Mean','Median','StdDev','Min','Max'}, 'RowNames', rowNames));

%% (b) Poisson regression
EasternTurkey = FMD(:,1);
lnCattle      = log(FMD(:,3));
X = [lnCattle, EasternTurkey];
y = FMD(:,2);

[beta, dev, stats] = glmfit(X, y, 'poisson', 'Link', 'log');

alpha = 0.05;
z  = norminv(1 - alpha/2);
lb = beta - z * stats.se;
ub = beta + z * stats.se;

fprintf('\n                 Poisson Regression Model\n');
rowlabel = {'Intercept', 'ln(Cattle)', 'EasternTurkey'};
fprintf('%-14s %8s %8s %8s %10s %10s\n', 'Variable', 'Coef', 'SE', 't-stat', 'CI Lower', 'CI Upper');
for i = 1:3
    fprintf('%-14s %8.3f %8.3f %8.3f %10.3f %10.3f\n', ...
        rowlabel{i}, beta(i), stats.se(i), stats.t(i), lb(i), ub(i));
end

%% (d)-(g) Expected counts and Pr(Y <= 2)
cases = { 'Eastern Turkey', 170080, 1;
          'Western Turkey', 157850, 0 };
for c = 1:size(cases,1)
    lambda = exp(beta(1) + beta(2)*log(cases{c,2}) + beta(3)*cases{c,3});
    fprintf('\n%s, Cattle = %d:\n', cases{c,1}, cases{c,2});
    fprintf('  Expected FMD count  lambda = %.4f\n', lambda);
    fprintf('  Pr(Y <= 2)                 = %.4f\n', poisscdf(2, lambda));
end
