% This directory contains test data and Matlab scripts to carry out an
% hemoglobin oxygen dissociation curve (ODC) using the graphic and/or
% kinetic method.

% 1. test data:  O2_pc_hb.csv
% 2. graphic method: HB_ODC_GRAPHIC_METHOD_O2_pc_hb.m
% 3. kinetic method: HB_ODC_KINETIC_METHOD_O2_pc_hb.m

% Additional required functions:
% reaction_0.m
% reaction.m
% simul_0.m
% simul_1.m
% simul_2.m
% simul_3.m

% To carry out the full kinetic method:
run HB_ODC_GRAPHIC_METHOD_O2_pc_hb.m 
run HB_ODC_KINETIC_METHOD_O2_pc_hb.m

% Execution of these scripts requires the Matlab Optimizazion Toolbox.