function [rI, cI, rU, cU] = imirawcall(data,sr,overlap)
%IMIRAWCALL   Calls IMIRAW.
%   IMIRAWCALL is a caller function for IMIRAW calculating raw entropy.
%   It plots and saves IMIRAW output.
%
%   [RI CI RU CU] = IMIRAWCALL(DATA,SR,OLP) calls IMIRAW with input
%   parameters DATA, SR (sampling rate) and OLP (overlap) and returns
%   IMIRAW output RI (real mutual information), CI (control mutual 
%   information), RU (real uncertainty coefficients) and CU (control
%   uncertainty coefficients).
%
%   See also IMIRAW.

% Directories
global DATAPATH
eg = '40';
resdir = [DATAPATH 'Ulbert\OITI_39_EEG_' eg '\MImap\'];
mm = pwd;

% Main
chno = size(data,2);    % prefiltering data is recommended (see IFILTER)!
lene = floor(size(data,1)/sr);
rI = zeros(chno,chno,(lene-1)*overlap+1);
rU = zeros(chno,chno,(lene-1)*overlap+1);
cI = zeros(chno,chno,(lene-1)*overlap+1);
cU = zeros(chno,chno,(lene-1)*overlap+1);
for x = 1:chno
    for y = x+1:chno
        disp([num2str(x) ' ' num2str(y)])
        ch1 = data(:,x)';
        ch2 = data(:,y)';
        [rHx,rHy,rHxy,rI(x,y,:),rU(x,y,:),rU(y,x,:),cHx,cHy,cHxy,cI(x,y,:),...
            cU(x,y,:),cU(y,x,:)] = imiraw(ch1,ch2,sr,overlap);
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
save(['MIrawmaps_EEG' eg '_filt01_40_overlap' num2str(overlap)],'rI','cI','rU','cU')
cd(mm)