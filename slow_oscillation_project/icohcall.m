function [nW rW, cW] = icohcall(data,sr)
%ICOHCALL   Calls ICOH.
%   ICOHCALL is a caller function for ICOH calculating wavelet coherence.
%   It plots and saves ICOH output.
%
%   See also ICOH.

% Directories
global DATAPATH
resdir = [DATAPATH 'Ulbert\WCoh\'];
mm = pwd;

% Main
chno = size(data,2);
lene = floor(size(data,1)/sr);
nW = zeros(chno,chno,lene);
rW = zeros(chno,chno,lene);
cW = zeros(chno,chno,lene);
for x = 1:chno
    for y = x+1:chno
        disp([num2str(x) ' ' num2str(y)])
        ch1 = data(:,x)';
        ch2 = data(:,y)';
        [rWCo cWCo] = icoh(ch1,ch2,sr);
        nWCo = rWCo - cWCo;
        srw1 = size(nWCo,1);
        srw2 = size(nWCo,2);
        nrWCo = reshape(nWCo,srw1*sr,srw2/sr);
        rrWCo = reshape(rWCo,srw1*sr,srw2/sr);
        crWCo = reshape(cWCo,srw1*sr,srw2/sr);
        nW(y,x,:) = mean(nrWCo);
        nW(x,y,:) = mean(nrWCo);
        rW(y,x,:) = mean(rrWCo);
        rW(x,y,:) = mean(rrWCo);
        cW(y,x,:) = mean(crWCo);
        cW(x,y,:) = mean(crWCo);
    end
end

% Plot
figure;
imagesc(mean(nW,3))
figure;
imagesc(mean(rW,3))
figure;
imagesc(mean(cW,3))

% Save
cd(resdir)
save WCOHmaps nW rW cW
cd(mm)