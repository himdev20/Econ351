%% Q3_probit.m
% ECO351 - Homework 1, Question 3
% Binary probit for support of marijuana legalization.
%   (a) Pr(y=1|x) = Phi(x'b),  log L = sum[ y log Phi(x'b) + (1-y) log(1-Phi(x'b)) ]
%   (b) MLE of the probit with covariates (1, x2, ..., x16)
%   (c) average covariate effect of age + 5 years
%   (d) average covariate effect of being a parent
%
% Data: data/marijuana.xlsx  (see data/README.md)

clear; clc; close all;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');

data = readtable(fullfile(dataDir, 'marijuana.xlsx'));
data.Properties.VariableNames = matlab.lang.makeValidName(strtrim(data.Properties.VariableNames));

% --- Recode (same scheme as Q2) ---------------------------------------
y = nan(height(data),1);
y(strcmpi(data.q85, 'Yes, legal'))  = 1;
y(strcmpi(data.q85, 'No, illegal')) = 0;

x2 = data.age;  x3 = data.hh1;
x4 = nan(height(data),1);
x4(strcmpi(data.pastUse,'Yes')) = 1;  x4(strcmpi(data.pastUse,'No')) = 0;
x5 = nan(height(data),1);
x5(strcmpi(data.sex,'Male')) = 1;     x5(strcmpi(data.sex,'Female')) = 0;
x6 = nan(height(data),1);
x6(strcmpi(data.parent,'Yes')) = 1;   x6(strcmpi(data.parent,'No')) = 0;
x7  = strcmpi(data.maritalStatus, 'Never been married');
x8  = ismember(data.maritalStatus, {'Divorced','Separated','Widowed'});
x9  = ismember(data.income, {'Less than 10000','10 to under 20000','20 to under 30000', ...
                             '30 to under 40000','40 to under 50000'});
x10 = ismember(data.income, {'50 to under 75000','75 to under 100000'});
x11 = ismember(data.educ, {'Less than HS','HS Incomplete','HS'});
x12 = ismember(data.educ, {'Some college','Associate Degree'});
x13 = strcmpi(data.race, 'White');
x14 = strcmpi(data.race, 'Black');
x15 = strcmpi(data.party, 'Democrat');
x16 = strcmpi(data.party, 'Republican');

X = [ones(height(data),1), x2, x3, x4, x5, x6, double(x7), double(x8), double(x9), ...
     double(x10), double(x11), double(x12), double(x13), double(x14), double(x15), double(x16)];
complete = all(~isnan(X), 2) & ~isnan(y);
Xc = X(complete,:);  yc = y(complete);

disp('================================================================');
disp('===== Question 3: Binary Probit (a-d) =====');
disp('================================================================');

%% (a) Probability expressions & likelihood
fprintf('Pr(y_i=1|X_i) = Phi(X_i*beta),  Pr(y_i=0|X_i) = 1 - Phi(X_i*beta)\n');
fprintf('l(beta) = sum[ y_i*log(Phi(X_i*beta)) + (1-y_i)*log(1-Phi(X_i*beta)) ]\n\n');

%% (b) Estimate the binary probit by MLE
[beta_hat, ~, stats] = glmfit(Xc(:,2:end), yc, 'binomial', 'link', 'probit');
se_hat = stats.se;  tvals = beta_hat ./ se_hat;

vars = {'Intercept','x2_age','x3_hh1','x4_pastUse','x5_sex_male','x6_parent', ...
        'x7_single','x8_postmarried','x9_poor','x10_middle','x11_HSandBelow', ...
        'x12_lessThanBachelors','x13_white','x14_black','x15_democrat','x16_republican'};
Tcoef = table(vars', round(beta_hat,3), round(se_hat,3), round(tvals,3), ...
    'VariableNames', {'term','coef','se','t'});
disp(Tcoef);

eta   = Xc * beta_hat;
p_hat = normcdf(eta);
ll    = sum(yc .* log(p_hat + 1e-12) + (1 - yc) .* log(1 - p_hat + 1e-12));
fprintf('Log-likelihood at MLE: %.6f\n', ll);

%% (c) Covariate effect: age + 5 years
APE_age5 = mean(normcdf(eta + 5*beta_hat(2)) - normcdf(eta));
fprintf('APE (Age + 5): %.4f\n', round(APE_age5, 4));

%% (d) Covariate effect: being a parent
X_parent1 = Xc; X_parent1(:,6) = 1;
X_parent0 = Xc; X_parent0(:,6) = 0;
ADE_parent = mean(normcdf(X_parent1*beta_hat) - normcdf(X_parent0*beta_hat));
fprintf('ADE (Parent): %.4f\n', round(ADE_parent, 4));
