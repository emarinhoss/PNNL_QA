% Early study: combine the flux from postprocessing.m (recon_004_MR_*) and
% restart_postprocessing.m (recon_MR_*) and plot the pooled mean with one
% standard deviation error bars.

clear all; close all; clc;

postprocessing
restart_postprocessing

tmean = zeros(1,2*frames+1);
tvar = tmean;

for k=1:10
    tmean = tmean + flux(k,:);
    tvar = tvar + flux(k,:).*flux(k,:);
end

for k=1:9
    tmean = tmean + flux2(k,:);
    tvar = tvar + flux2(k,:).*flux2(k,:);
end

tmean = tmean /19;

tvar = tvar/19 - tmean.*tmean;

errorbar(0:2*frames,tmean,sqrt(tvar),'r')

