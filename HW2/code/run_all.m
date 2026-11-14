%% run_all.m -- run every Homework 2 script
% Questions 1 and 4 are analytical; see report/HW2-solutions.pdf.
here = fileparts(mfilename('fullpath'));
for s = {'Q2_cola_mnl', 'Q3_fmd_poisson'}
    fprintf('\n\n######## Running %s ########\n', s{1});
    run(fullfile(here, [s{1} '.m']));
end
