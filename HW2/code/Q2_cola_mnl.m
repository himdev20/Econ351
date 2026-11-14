%% Q2_cola_mnl.m
% ECO351 - Homework 2, Question 2
% Multinomial logit for soft-drink brand choice (Pepsi, 7-Up, Coke Classic).
% Data are stacked: 3 rows per customer (1 = Pepsi, 2 = 7-Up, 3 = Coke).
%   (b) fit MNL of choice on price (Coke = base category)
%   (c) predicted probabilities at Pepsi = $1.00, 7-Up = $1.25, Coke = $1.10
%   (d) change in probabilities from a $0.10 price increase for each brand
%
% Data: data/cola.xlsx  (see data/README.md)

clc; clear;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
data    = readtable(fullfile(dataDir, 'cola.xlsx'));

N = length(unique(data.id));
J = 3;

% Reshape stacked data to one row per customer
choiceMat = reshape(data.choice, [J, N])';
[~, y]    = max(choiceMat, [], 2);              % y in {1,2,3}
priceMat  = reshape(data.price,  [J, N])';

% Prices relative to the base alternative (Coke)
X = [priceMat(:,1) - priceMat(:,3), priceMat(:,2) - priceMat(:,3)];

%% (b) Estimate
[B, dev, stats] = mnrfit(X, y, 'model', 'nominal');

disp('Coefficient estimates:');  disp(B);
disp('Standard errors:');        disp(stats.se);
disp('t-statistics:');           disp(B ./ stats.se);

%% (c) Predicted probabilities at the given prices
priceCoke = 1.10;  pricePepsi = 1.00;  price7Up = 1.25;
Xnew = [pricePepsi - priceCoke, price7Up - priceCoke];   % Coke is the base, so its diff is 0
prob = mnrval(B, Xnew, 'model', 'nominal');

disp('Predicted choice probabilities:');
fprintf('P(Pepsi) = %.4f\nP(7-Up)  = %.4f\nP(Coke)  = %.4f\n', prob(1), prob(2), prob(3));

%% (d) Effect of a 10-cent price increase for each brand
scenarios = { 'Pepsi', [1.10, 1.25, 1.10];      % [pricePepsi, price7Up, priceCoke]
              '7-Up',  [1.00, 1.35, 1.10];
              'Coke',  [1.00, 1.25, 1.20] };

for s = 1:size(scenarios,1)
    p    = scenarios{s,2};
    Xs   = [p(1) - p(3), p(2) - p(3)];
    prob_s = mnrval(B, Xs, 'model', 'nominal');
    fprintf('\n-------- Increasing price of %s by 10 cents --------\n', scenarios{s,1});
    fprintf('New probabilities  [Pepsi 7-Up Coke]: %s\n', mat2str(prob_s, 4));
    fprintf('Change in prob.    [Pepsi 7-Up Coke]: %s\n', mat2str(prob_s - prob, 4));
end
