% Postprocessing
% Calculates the mean and variance in the reconnection
% flux for the different cases.

clear all; clc;

d = dir('recon_004_MR_*');
numdir = size(d);

frames = 40;

flux = zeros(numdir(1),frames+1);
fmean= zeros(1,frames+1);
fvar = fmean;

for k=1:numdir(1)
    folder = d(k).name
    for l=0:frames
        [output, x, y] = load_data_new([folder '/ssrecon_wv'],'qnew',l);
        mid=ceil(0.5*length(y));
        flux(k,l+1) = sum(abs(output(mid,:,15)));
    end
end

dx = abs(x(2)-x(1));
Ly = 12.8;
flux = 0.2*flux/(2*Ly);

for k=1:numdir(1)
    flux(k,:) = flux(k,:)/flux(k,1);
    fmean = fmean + flux(k,:);
    fvar  = fvar + flux(k,:).*flux(k,:);
end

fmean = fmean / numdir(1);
fvar = fvar / numdir(1) - fmean.*fmean;

errorbar(0:frames,fmean,fvar,'b'), hold on