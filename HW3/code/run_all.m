%% run_all.m -- run every Homework 3 script
% Question 1 is analytical; see report/HW3-solutions.pdf.
here = fileparts(mfilename('fullpath'));
for s = {'Q2_tobit', 'Q3_heckman', 'Q4_treatment_effect'}
    fprintf('\n\n######## Running %s ########\n', s{1});
    run(fullfile(here, [s{1} '.m']));
end
