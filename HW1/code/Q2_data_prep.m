%% Q2_data_prep.m
% ECO351 - Homework 1, Question 2
% Pew Research Center (March 2013) survey on marijuana legalization, n = 1240.
% Recodes the raw survey responses into the response y and covariates x2..x16
% used by the probit (Q3) and logit (Q4) models, and prints frequency tables.
%
% Output: HW1/output/prepared_data.mat and prepared_data.csv
% Data:   data/marijuana.xlsx  (see data/README.md)

clear; clc;

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', '..', 'data');
outDir  = fullfile(here, '..', 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

data = readtable(fullfile(dataDir, 'marijuana.xlsx'));
N = height(data);

%% (a) Response variable: y = 1 if 'Yes, legal', 0 if 'No, illegal'
y = nan(N,1);
y(strcmpi(data.q85, 'Yes, legal'))  = 1;
y(strcmpi(data.q85, 'No, illegal')) = 0;
disp('--- (a) Response variable y ---');
disp(array2table(tabulate(y), 'VariableNames', {'Category','Count','Percent'}));

%% (b) Age (x2) and household size (x3)
x2 = data.age;
x3 = data.hh1;
disp('--- (b) Age and Household ---');
disp(table(mean(x2,'omitnan'), std(x2,'omitnan'), mean(x3,'omitnan'), std(x3,'omitnan'), ...
    'VariableNames', {'Mean_Age','Std_Age','Mean_HH1','Std_HH1'}));

%% (c) Past use (x4)
x4 = nan(N,1);
x4(strcmpi(data.pastUse, 'Yes')) = 1;
x4(strcmpi(data.pastUse, 'No'))  = 0;
disp('--- (c) Past Use ---');
disp(array2table(tabulate(x4), 'VariableNames', {'Category','Count','Percent'}));

%% (d) Sex (x5 = 1 if male) and parent (x6)
x5 = nan(N,1);
x5(strcmpi(data.sex, 'Male'))   = 1;
x5(strcmpi(data.sex, 'Female')) = 0;
x6 = nan(N,1);
x6(strcmpi(data.parent, 'Yes')) = 1;
x6(strcmpi(data.parent, 'No'))  = 0;
disp('--- (d) Sex ---');
disp(array2table(tabulate(x5), 'VariableNames', {'Category','Count','Percent'}));
disp('--- (d) Parent ---');
disp(array2table(tabulate(x6), 'VariableNames', {'Category','Count','Percent'}));

%% (e) Marital status: single (x7), post-married (x8), couple (base)
x7 = strcmpi(data.maritalStatus, 'Never been married');
x8 = ismember(data.maritalStatus, {'Divorced','Separated','Widowed'});
couple = ~(x7 | x8);
print_categories('(e) Marital Status', {'Single';'Post-Married';'Couple'}, [x7, x8, couple]);

%% (f) Income: poor (x9), middle (x10), rich (base)
poor   = ismember(data.income, {'Less than 10000','10 to under 20000','20 to under 30000', ...
                                '30 to under 40000','40 to under 50000'});
middle = ismember(data.income, {'50 to under 75000','75 to under 100000'});
rich   = ~(poor | middle);
x9 = poor; x10 = middle;
print_categories('(f) Income', {'Poor';'Middle';'Rich'}, [poor, middle, rich]);

%% (g) Education: HS and below (x11), less than bachelors (x12), bachelors+ (base)
HSandBelow        = ismember(data.educ, {'Less than HS','HS Incomplete','HS'});
lessThanBachelors = ismember(data.educ, {'Some college','Associate Degree'});
BachelorsAndAbove = ~(HSandBelow | lessThanBachelors);
x11 = HSandBelow; x12 = lessThanBachelors;
print_categories('(g) Education', {'HS and Below';'Less Than Bachelors';'Bachelors and Above'}, ...
    [HSandBelow, lessThanBachelors, BachelorsAndAbove]);

%% (h) Race: white (x13), black (x14), all other (base)
white    = strcmpi(data.race, 'White');
black    = strcmpi(data.race, 'Black');
allOther = ~(white | black);
x13 = white; x14 = black;
print_categories('(h) Race', {'White';'Black';'All Other'}, [white, black, allOther]);

%% (i) Party: democrat (x15), republican (x16), independent/others (base)
democrat          = strcmpi(data.party, 'Democrat');
republican        = strcmpi(data.party, 'Republican');
independentOthers = ~(democrat | republican);
x15 = democrat; x16 = republican;
print_categories('(i) Party', {'Democrat';'Republican';'Independent/Others'}, ...
    [democrat, republican, independentOthers]);

%% Save the prepared design matrix
X = double([x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16]);
varNames = strcat("x", string(2:16));
save(fullfile(outDir, 'prepared_data.mat'), 'y', 'X', 'varNames');
writetable(array2table([y X], 'VariableNames', ['y', varNames]), ...
           fullfile(outDir, 'prepared_data.csv'));
fprintf('\nPrepared dataset saved to %s\n', outDir);

%% ---------------------------------------------------------------------
function print_categories(titleStr, labels, D)
    counts = sum(D, 1);
    perc   = counts / size(D,1) * 100;
    disp(['--- ' titleStr ' ---']);
    disp(table(labels, counts', perc', 'VariableNames', {'Category','Count','Percent'}));
end
