%% run_all.m -- run every Homework 1 script in order
% Q4 depends on the prepared data written by Q2, so order matters.
here = fileparts(mfilename('fullpath'));
scripts = {'Q1_igg', 'Q2_data_prep', 'Q3_probit', 'Q4_logit', 'Q5_ordinal_probit'};
for s = scripts
    fprintf('\n\n######## Running %s ########\n', s{1});
    run(fullfile(here, [s{1} '.m']));
end
