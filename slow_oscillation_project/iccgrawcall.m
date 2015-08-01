function [rI, cI] = iccgrawcall(data,sr,overlap)
%ICCGRAWCALL   Calls ICCGRAW.
%   ICCGRAWCALL is a caller function for ICCGRAW calculating
%   cross-correlation with shuffled control. It plots and saves ICCGRAW
%   output.
%
%   [RI CI] = ICCGRAWCALL(DATA,SR,OLP) calls ICCGRAW with input parameters
%   DATA, SR (sampling rate) and OLP (overlap) and returns ICCGRAW output
%   RI (real cross-correlation) and CI (control cross-correlation).
%
%   See also ICCGRAW.

% Directories
global DATAPATH
eg = '298';
resdir = [DATAPATH 'Ulbert\OITI_N1_EEG_' eg '\CCGmap\'];
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
        [rI(x,y,:),cI(x,y,:)] = iccgraw(ch1,ch2,sr,overlap);
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
save(['MIrawmaps_EEG' eg '_filt01_40_overlap' num2str(overlap)],'rI','cI')
cd(mm)