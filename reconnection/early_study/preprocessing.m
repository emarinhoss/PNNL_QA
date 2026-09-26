% Early study: Monte Carlo setup for the reconnection runs.
%
% Creates `number` run folders recon_004_MR_<MI/ME>_c0_<c> with randomly
% sampled ion-to-electron mass ratio (uniform on [95, 105], ME = 0.01) and
% speed of light (uniform on [1, 3]), writes ssrecon_wv.pin from input.py,
% preprocesses it with $wxpp and moves ssrecon_wv.inp into the folder.

clear all; close all; clc;

% desired number of runs to be executed
number = 300;

% Run parameter that is similar in all runs
sss = fopen('input.py','r');
info = fscanf(sss,'%c');
fclose(sss);

% Generate ramdom data for Mass ratio and speed of light
for k=1:number
    % uniform varying values of mass ration from
    % a to b given a certain value of specied value
    % electron mass
    a = 95; b = 105;
    v = a + (b-a).*rand(1);
    me = 0.01;
    mi = v*me;
    
    % uniform varying values of the speed of light
    d = 1; e = 3;
    c = d + (e-d).*rand(1);
    light = ['LIGHT_SPEED = ' num2str(c)];
    masse = ['ME = ' num2str(me)];
    massi = ['MI = ' num2str(mi)];
    
    % write new input file
    out = fopen('ssrecon_wv.pin','w');
    fprintf(out,'# -*- python -*- \n');
    fprintf(out,'# The following parameters has been randomly generated. \n');
    fprintf(out,light);
    fprintf(out,'\n');
    fprintf(out,masse);
    fprintf(out,'\n');
    fprintf(out,massi);
    fprintf(out,'\n');
    fprintf(out,'# -- End of randomly generated data. -- \n');
    fprintf(out,'\n');
    fprintf(out,info);
    fclose(out);
    
    % create a folder for each run
    folder= ['recon_004_MR_' num2str(v) '_c0_' num2str(c)];
    mkdir(folder)
    
    % parse input file for run
    system('$wxpp -i ssrecon_wv.pin');
    
    % muve input file into run folder
    system(['mv ssrecon_wv.inp ' folder]);
end