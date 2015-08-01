function b_3kwentropy_control2
%4KWENTROPY_CONTROL2   Wavelet entropy calculation.
%   4KWENTROPY_CONTROL2 calculates wavelet magnitude entropies for 2MEGAWAVERUN_ENTRCTRL output
%   files. It calls ENTR3 for the raw entropy calculation. 3KWENTROPY_CONTROL2 computes entropy
%   for eeg wavelet and delayed unit wavelet (10 seconds delay). Outputs:
%       Rhxabs: unit wavelet magnitude entropy
%       Rhyabs: eeg wavelet magnitude entropy
%       Rhxyabs: combined entropy
%       Rixyabs: mutual information
%       Rixynormabs: normalized mutual information
%       Rrelshanyabs: relative Shannon entropy for eeg
%       Rrelshanxabs: relative Shannon entropy for unit
%       Rhxcyabs: conditional entropy (H(unit|eeg))
%       Rhycxabs: conditional entropy (H(eeg|unit))
%       Ruxyabs: uncertainity coefficient (eeg->unit)
%       Ruyxabs: uncertainity coefficient (unit->eeg)
%
%   See also ENTR3 and 2MEGAWAVERUN_ENTRCTRL.

% Input argument check
error(nargchk(0,0,nargin))

% Define directories
global DATADIR
global DATAPATH
wh_from = [DATAPATH 'Entropy2\control\megawave\'];
wh_to = [DATAPATH 'Entropy4\control\delayed_unit_10sec\'];
files = dir(wh_from);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
    end
end
files2 = files2(2:end);
sf = length(files2);

% Progress indicator
wb = waitbar(0,'Running 4KWENTROPY CONTROL2...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
for i = 1:11    % CELL CYCLE
    fnm = files2(i).name;
    ffnm = [wh_from '\' fnm];
    load(ffnm);
    
    load([wh_to 'settings\ScaleVector']);   % load scale vector
    fnd = find(f>6);    % frequency band bounderies
    pwind1 = fnd(end);
    fnd = find(f<2.5);
    pwind2 = fnd(1);
    
    abs1 = wavea(pwind1:pwind2,:);  % eeg
    abs2 = wavef(pwind1:pwind2,:);  % delayed unit (10 sec.)
    clear wavea waveb wavec waved wavee wavef
    [k1 k2] = size(abs1);
    
    winlen = 5000;   % window size
    entloopend = k2 / winlen;
    
    Rhxabs = []; Rhyabs = []; Rhxyabs = [];     % preallocating output variables
    Rhycxabs = []; Rhxcyabs = [];
    Rixyabs = []; Ruxyabs = []; Ruyxabs = [];
    Rixynormabs = [];
    Rrelshanyabs = []; Rrelshanxabs = [];
    
% ABS LOOP
    for entloop = 1:entloopend*10-10        % WINDOW LOOP
        inx1 = (entloop - 1) * winlen / 10 + 1;  % Note: overlaping windows!
        inx2 = inx1 + winlen - 1;
        
        y1 = abs1(:,inx1:inx2);   % eeg wavelet "cells"
        y2 = abs2(:,inx1:inx2);   % unit wavelet "cells"
        leny = numel(y1);
        hb = fix(exp(0.626+0.4*log(leny-1)));   % number of bins
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entr3(abs1,abs2,y1,hb,y2,hb); % MAIN
        Ixynorm = Ixy / log(leny);
        relshany = Hy / log(hb);
        relshanx = Hx / log(hb);
        
        Rhxabs = [Rhxabs Hx]; 
        Rhyabs = [Rhyabs Hy]; 
        Rhxyabs = [Rhxyabs Hxy]; 
        Rhycxabs = [Rhycxabs Hycx]; 
        Rhxcyabs = [Rhxcyabs Hxcy]; 
        Rixyabs = [Rixyabs Ixy]; 
        Ruxyabs = [Ruxyabs Uxy]; 
        Ruyxabs = [Ruyxabs Uyx]; 
        Rixynormabs = [Rixynormabs Ixynorm]; 
        Rrelshanyabs = [Rrelshanyabs relshany];
        Rrelshanxabs = [Rrelshanxabs relshanx];
    end     % end of abs (window) loop
    
% Save
    fs = findstr(fnm,'_');
    fn = fnm(fs(2)+1:end-4);
    fnts = [wh_to fn '_10SEC'];
    save(fnts,'Rhxabs','Rhyabs','Rhxyabs','Rixyabs','Rixynormabs','Rrelshanyabs','Rrelshanxabs',...
        'Rhxcyabs','Rhycxabs','Ruxyabs','Ruyxabs');
    clear abs1 abs2
    
    waitbar(i/sf)   % progress indicator
end     % end of cell cycle
close(wb)