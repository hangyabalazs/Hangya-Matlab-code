function [rI, cI, rU, cU] = ientrycall(data,sr)
%IENTRYCALL   Calls IENTRY.
%   IENTRYCALL is a caller function for IENTRY calculating wavelet entropy.
%   It plots and saves IENTRY output.
%
%   See also IENTRY.

% Directories
global DATAPATH
resdir = [DATAPATH 'Ulbert\OITI_N1_EEG_148a\MImap\'];
mm = pwd;

% Main
chno = size(data,2);
lene = floor(size(data,1)/sr);
rI = zeros(chno,chno,lene);
rU = zeros(chno,chno,lene);
cI = zeros(chno,chno,lene);
cU = zeros(chno,chno,lene);
for x = 1:chno
    for y = x+1:chno
        disp([num2str(x) ' ' num2str(y)])
        ch1 = data(:,x)';
        ch2 = data(:,y)';
        [rHx,rHy,rHxy,rI(x,y,:),rU(x,y,:),rU(y,x,:),cHx,cHy,cHxy,cI(x,y,:),...
            cU(x,y,:),cU(y,x,:)] = ientry(ch1,ch2,sr);
        rI(y,x,:) = rI(x,y,:);
        cI(y,x,:) = cI(x,y,:);
    end
end

% Plot
figure;
imagesc(mean(rI,3))
figure;
imagesc(mean(cI,3))

figure;
imagesc(mean(rU,3))
figure;
imagesc(mean(cU,3))

% Save
cd(resdir)
save MImaps_EEG148 rI cI rU cU
cd(mm)