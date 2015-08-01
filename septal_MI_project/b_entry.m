function b_entry
%ENTRY   Wavelet entropy calculation.
%   ENTRY calculates wavelet magnitude entropies for MEGAWAVERUN_ENTRCTRL2 output files. 
%   It calls ENTRY_SUB for the raw entropy calculation. Being a third generation entropy function,
%   it computes entropy for the theta band only and normalizes for the whole theta band. It
%   uses 5 sec. long, 80 % overlaping windows. Outputs:
%       aHx: unit wavelet magnitude entropy
%       aHy: eeg wavelet magnitude entropy
%       aHxy: combined entropy
%       aIxy: mutual information
%       aIxynorm: normalized mutual information
%       aRelHy: relative Shannon entropy for eeg
%       aRelHx: relative Shannon entropy for unit
%       aHxcy: conditional entropy (H(unit|eeg))
%       aHycx: conditional entropy (H(eeg|unit))
%       aUxy: uncertainity coefficient (eeg->unit)
%       aUyx: uncertainity coefficient (unit->eeg)
%
%   See also ENTRY_SUB and MEGAWAVERUN_ENTRCTRL2.

% Input argument check
error(nargchk(0,0,nargin))

% Define directories
global DATADIR
global DATAPATH
wh_from = [DATAPATH 'Entropy2\control\megawave\'];
wh_to = [DATAPATH 'Entropy3\temp2\'];
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
wb = waitbar(0,'Running ENTRY...','Position',[360 250 275 50]);    %Progress indicator
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
    abs2 = waveb(pwind1:pwind2,:);  % unit
    clear wavea waveb wavec waved wavee wavef
    [k1 k2] = size(abs1);
    
    winlen = 5000;   % window size
    maxi = k2 / winlen;
    
    aHx = []; aHy = []; aHxy = [];     % preallocating output variables
    aHycx = []; aHxcy = [];
    aIxy = []; aUxy = []; aUyx = [];
    aIxynorm = [];
    aRelHy = []; aRelHx = [];
    
% ABS LOOP
    for i = 1:maxi*10-10        % WINDOW LOOP
        inx1 = (i - 1) * winlen / 10 + 1;  % Note: overlaping windows!
        inx2 = inx1 + winlen - 1;
        
        y1 = abs1(:,inx1:inx2);   % eeg wavelet "cells"
        y2 = abs2(:,inx1:inx2);   % unit wavelet "cells"
        numy = numel(y1);
        bno = fix(exp(0.626+0.4*log(numy-1)));   % number of bins
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entry_sub(abs1,abs2,y1,bno,y2,bno); % MAIN
        Ixynorm = Ixy / log(numy);
        RelHy = Hy / log(bno);      % relative Shannon entropy
        RelHx = Hx / log(bno);
        
        aHx = [aHx Hx]; 
        aHy = [aHy Hy]; 
        aHxy = [aHxy Hxy]; 
        aHycx = [aHycx Hycx]; 
        aHxcy = [aHxcy Hxcy]; 
        aIxy = [aIxy Ixy]; 
        aUxy = [aUxy Uxy];
        aUyx = [aUyx Uyx]; 
        aIxynorm = [aIxynorm Ixynorm]; 
        aRelHy = [aRelHy RelHy];
        aRelHx = [aRelHx RelHx];
    end     % end of abs (window) loop
    
% Save
    fs = findstr(fnm,'_');
    fn = fnm(fs(2)+1:end-4);
    fnts = [wh_to fn '_ENTROPY'];
    save(fnts,'aHx','aHy','aHxy','aIxy','aIxynorm','aRelHy','aRelHx',...
        'aHxcy','aHycx','aUxy','aUyx');
    clear abs1 abs2
    
    waitbar(i/sf)   % progress indicator
end     % end of cell cycle
close(wb)