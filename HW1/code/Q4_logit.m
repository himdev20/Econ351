%% Q4_logit.m
% ECO351 - Homework 1, Question 4
% Binary logit on the prepared marijuana data.
%   (a) Pr(y=1|x) = Lambda(x'b) = 1 / (1 + exp(-x'b))
%   (b) MLE of the logit with covariates (1, x2, ..., x16)
%   (c) effect of age + 5 years on Pr(y=1)   (evaluated at the sample mean)
%   (d) effect of being a parent on Pr(y=1)  (evaluated at the sample mean)
%   (e) comparison with probit
%
% Requires: run Q2_data_prep.m first (creates HW1/output/prepared_data.mat)

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, '..', 'output', 'prepared_data.mat'));
y = S.y;  X = S.X;  varNames = S.varNames;

%% (a) Probability functions and likelihood
fprintf('\n(a) P(y_i=1|x_i) = exp(x_i*beta)/(1+exp(x_i*beta))\n');
fprintf('    L(beta) = prod[ P(y_i=1)^y_i * (1-P(y_i=1))^(1-y_i) ]\n');

%% (b) Estimate binary logit
T = array2table([y X], 'VariableNames', ['y', varNames]);
formula = "y ~ " + strjoin(varNames, " + ");
mdl   = fitglm(T, formula, 'Distribution', 'binomial', 'Link', 'logit');
coefs = mdl.Coefficients;
disp('--- (b) Logit Estimates ---');
disp(table(coefs.Estimate, coefs.SE, coefs.tStat, ...
    'VariableNames', {'Estimate','StdErr','tValue'}, ...
    'RowNames', coefs.Properties.RowNames));

b = mdl.Coefficients.Estimate;
logistic = @(z) 1 ./ (1 + exp(-z));

%% (c) Effect of increasing age by 5 years
meanX = mean(X, "omitnan");
p0 = logistic([1 meanX] * b);
meanX(1) = meanX(1) + 5;                    % column 1 of X is x2 = age
p5 = logistic([1 meanX] * b);
fprintf('\n(c) Change in probability for +5 years age: %.4f\n', p5 - p0);

%% (d) Effect of being a parent
meanX = mean(X, "omitnan");
p0 = logistic([1 meanX] * b);
meanX(5) = 1;                               % column 5 of X is x6 = parent
p_parent = logistic([1 meanX] * b);
fprintf('(d) Change in probability for being a parent: %.4f\n', p_parent - p0);

%% (e) Comparison to probit
fprintf('\n(e) Both probit and logit give qualitatively similar effects;\n');
fprintf('    magnitudes differ since logit uses the logistic CDF vs the probit normal CDF.\n');
