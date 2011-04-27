clc;

d = dir('recon_MR_*');
numdir = size(d);

frames = 40;

flux2 = zeros(numdir(1),2*frames+1);
fmean2= zeros(1,2*frames+1);
fvar2 = fmean2;

for k=1:numdir(1)
    folder = d(k).name;
    for l=0:frames
        [output1, x, y] = load_data_new([folder '/ssrecon_wv'],'qnew',l);
        [output2, x, y] = load_data_new([folder '/ssrecon_wv_r'],'qnew',l+frames);
        mid=ceil(0.5*length(y));
        flux2(k,l+1) = sum(abs(output1(mid,:,15)));
        flux2(k,l+1+frames) = sum(abs(output2(mid,:,15)));
    end
end

dx = abs(x(2)-x(1));
Ly = 12.8;
flux2 = 0.2*flux2/(2*Ly);

for k=1:numdir(1)
    flux2(k,:) = flux2(k,:)/flux2(k,1);
    fmean2 = fmean2 + flux2(k,:);
    fvar2  = fvar2 + flux2(k,:).*flux2(k,:);
end

fmean2 = fmean2 / numdir(1);
fvar2 = fvar2 / numdir(1) - fmean2.*fmean2;

errorbar(0:2*frames,fmean2,fvar2,'g')