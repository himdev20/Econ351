%% Q3_murder_panel.m
% ECO351 - Homework 4, Question 3
% US state panel: murder rates and executions (years 87, 90, 93).
%   mrdrte_it = eta_t + b1*exec_it + b2*unemp_it + c_i + e_it
%   (b) pooled OLS, 1990 & 1993            (c) first differences, 1990 -> 1993
%   (d) robust vs non-robust SEs           (e) state with most executions in 1993
%   (f) first differences without Texas    (g) within (FE) estimator, all years,
%                                              SEs clustered by state
%
% Data: data/Murder.xlsx  (sheet "Data": state, year, mrdrte, exec, unem)

clear; close all; clc;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
TBL = readtable(fullfile(dataDir, 'Murder.xlsx'), 'Sheet', 'Data');

state  = string(TBL.state);
year   = double(TBL.year);
mrdrte = double(TBL.mrdrte);
exec   = double(TBL.exec);
unemp  = double(TBL.unem);
data   = table(state, year, mrdrte, exec, unemp);

fprintf('Loaded %d observations across %d states and %d years: %s\n', ...
    height(data), numel(unique(state)), numel(unique(year)), mat2str(unique(year)'));

%% (b) Pooled OLS using 1990 and 1993
df23 = data(ismember(data.year, [90, 93]), :);
D93  = double(df23.year == 93);
X_pooled = [ones(height(df23),1), D93, df23.exec, df23.unemp];
res_pooled = ols_infer(df23.mrdrte, X_pooled);

fprintf('\n(b) Pooled OLS (1990 & 1993):\n');
names = {'Intercept','Year1993','exec','unemp'};
for j = 1:numel(names)
    fprintf('%10s: %9.4f  (robust SE = %6.4f)  p = %6.4f\n', ...
        names{j}, res_pooled.beta(j), res_pooled.se_HC1(j), res_pooled.p(j));
end

%% (c) Fixed effects via first differencing (1993 - 1990)
states = unique(df23.state);
diffs  = table('Size', [0 4], 'VariableTypes', {'string','double','double','double'}, ...
               'VariableNames', {'state','d_mrdrte','d_exec','d_unemp'});
for i = 1:numel(states)
    s  = states(i);
    r0 = df23(df23.state == s & df23.year == 90, :);
    r1 = df23(df23.state == s & df23.year == 93, :);
    if ~isempty(r0) && ~isempty(r1)
        diffs = [diffs; {s, r1.mrdrte - r0.mrdrte, r1.exec - r0.exec, r1.unemp - r0.unemp}]; %#ok<AGROW>
    end
end
res_fd = ols_infer(diffs.d_mrdrte, [diffs.d_exec, diffs.d_unemp]);     % no intercept

fprintf('\n(c) First difference (1993 - 1990):\n');
names_fd = {'d_exec','d_unemp'};
for j = 1:2
    fprintf('%10s: %9.4f  (robust SE = %6.4f)  p = %6.4f\n', ...
        names_fd{j}, res_fd.beta(j), res_fd.se_HC1(j), res_fd.p(j));
end

%% (d) Non-robust vs robust SE (pooled OLS, exec)
fprintf('\n(d) Pooled OLS, exec: non-robust SE = %.4f ; robust HC1 SE = %.4f\n', ...
    res_pooled.se_ols(3), res_pooled.se_HC1(3));

%% (e) State with the most executions in 1993
df93 = sortrows(data(data.year == 93, :), 'exec', 'descend');
fprintf('\n(e) Top state:    %s | executions = %.0f\n', df93.state(1), df93.exec(1));
fprintf('    Second state: %s | executions = %.0f\n', df93.state(2), df93.exec(2));
fprintf('    Difference:   %.0f executions\n', df93.exec(1) - df93.exec(2));

%% (f) First differences excluding Texas
diffs_noTX  = diffs(diffs.state ~= "TX", :);
res_fd_noTX = ols_infer(diffs_noTX.d_mrdrte, [diffs_noTX.d_exec, diffs_noTX.d_unemp]);

fprintf('\n(f) First difference, excluding Texas:\n');
for j = 1:2
    fprintf('%10s: beta = %9.4f | non-robust SE = %7.4f | robust SE = %7.4f | p = %6.4f\n', ...
        names_fd{j}, res_fd_noTX.beta(j), res_fd_noTX.se_ols(j), res_fd_noTX.se_HC1(j), res_fd_noTX.p(j));
end

%% (g) Within (FE) estimator, all years, cluster-robust SEs by state
N  = height(data);
uy = unique(data.year);
yearDummies = zeros(N, numel(uy) - 1);           % base year = first year
for j = 2:numel(uy)
    yearDummies(:, j-1) = (data.year == uy(j));
end
Xraw = [data.exec, data.unemp, yearDummies];
Yall = data.mrdrte;

groups = unique(data.state);
G = numel(groups);
Y_within  = zeros(N,1);
X_within  = zeros(N, size(Xraw,2));
clusterid = zeros(N,1);
for i = 1:G
    idx = (data.state == groups(i));
    Y_within(idx)   = Yall(idx)   - mean(Yall(idx));
    X_within(idx,:) = Xraw(idx,:) - mean(Xraw(idx,:), 1);
    clusterid(idx)  = i;
end

beta_within = (X_within' * X_within) \ (X_within' * Y_within);
se_cluster  = sqrt(diag(clusterV(Y_within, X_within, clusterid)));
zstat       = beta_within ./ se_cluster;
p_cluster   = 2 * (1 - normcdf(abs(zstat)));

fprintf('\n(g) Within estimator, all years (SE clustered by state, G = %d):\n', G);
names_within = [{'exec','unemp'}, arrayfun(@(y) sprintf('year%d', y), uy(2:end)', 'UniformOutput', false)];
for j = 1:numel(beta_within)
    fprintf('%10s: %9.4f  (cluster SE = %6.4f)  p = %6.4f\n', ...
        names_within{j}, beta_within(j), se_cluster(j), p_cluster(j));
end

%% Summary: exec coefficient across methods
fprintf('\n--- Summary: exec coefficient across methods ---\n');
fprintf('Pooled OLS (1990 & 1993) : %8.4f  (SE = %.4f, p = %.4f)\n', res_pooled.beta(3),  res_pooled.se_HC1(3),  res_pooled.p(3));
fprintf('FD (1990 -> 1993)        : %8.4f  (SE = %.4f, p = %.4f)\n', res_fd.beta(1),      res_fd.se_HC1(1),      res_fd.p(1));
fprintf('FD without Texas         : %8.4f  (SE = %.4f, p = %.4f)\n', res_fd_noTX.beta(1), res_fd_noTX.se_HC1(1), res_fd_noTX.p(1));
fprintf('FE, all years            : %8.4f  (SE = %.4f, p = %.4f)\n', beta_within(1),      se_cluster(1),         p_cluster(1));

%% ---------------------------------------------------------------------
%% Local functions
%% ---------------------------------------------------------------------
function out = ols_infer(Y, X)
    % OLS with conventional and HC1 heteroskedasticity-robust standard errors
    [n, k] = size(X);
    XXinv = inv(X' * X);
    beta  = XXinv * (X' * Y);
    e     = Y - X * beta;
    s2    = (e' * e) / (n - k);
    Var_beta = s2 * XXinv;
    Var_HC1  = XXinv * (X' * diag(e.^2) * X) * XXinv * (n / (n - k));

    out.beta     = beta;
    out.se_ols   = sqrt(diag(Var_beta));
    out.se_HC1   = sqrt(max(0, diag(Var_HC1)));
    out.t        = beta ./ out.se_HC1;
    out.p        = 2 * (1 - tcdf(abs(out.t), n - k));
    out.resid    = e;
    out.Var_beta = Var_beta;
    out.Var_HC1  = Var_HC1;
end

function Vcl = clusterV(Y, X, clust)
    % Cluster-robust (sandwich) variance of the OLS estimator
    k = size(X,2);
    beta = (X' * X) \ (X' * Y);
    e = Y - X * beta;
    groups = unique(clust);
    XXinv = inv(X' * X);
    S = zeros(k);
    for g = 1:numel(groups)
        idx = (clust == groups(g));
        sg  = X(idx,:)' * e(idx);
        S   = S + sg * sg';
    end
    Vcl = XXinv * S * XXinv;
end
