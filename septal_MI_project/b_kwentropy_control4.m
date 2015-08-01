function b_kwentropy_control4
%KWENTROPY_CONTROL4   Wavelet entropy calculation.
%   KWENTROPY_CONTROL4 calculates wavelet magnitude entropies for MEGAWAVERUN_ENTRCTRL2 output
%   files. It calls ENTR2 for the raw entropy calculation. KWENTROPY_CONTROL4 computes entropy
%   for eeg wavelet and random data wavelet. Outputs:
%       Rhxabs: unit wavelet magnitude entropy
%       Rhyabs: eeg wavelet magnitude entropy
%       Rhxyabs: combined entropy
%       Rixyabs: mutual information
%       Rixynormabs: normalized mutual information
%       Rrelshanabs: relative Shannon entropy for eeg
%       Rhxcyabs: conditional entropy (H(unit|eeg))
%       Rhycxabs: conditional entropy (H(eeg|unit))
%       Ruxyabs: uncertainity coefficient (eeg->unit)
%       Ruyxabs: uncertainity coefficient (unit->eeg)
%
%   See also ENTR2 and MEGAWAVERUN_ENTRCTRL2.

% Input argument check
error(nargchk(0,0,nargin))

% Define directories
global DATADIR
global DATAPATH
wh_from = [DATAPATH 'Entropy\control\megawave\'];
wh_to = [DATAPATH 'Entropy\control\random_vs_eeg\'];
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
wb = waitbar(0,'Running KWENTROPY CONTROL4...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
for i = 1:10    % CELL CYCLE
    fnm = files2(i).name;
    ffnm = [wh_from '\' fnm];
    load(ffnm);
    
    abs1 = wavea;  % eeg
    clear wavea;
    abs2 = wavec;  % random
    clear wavec;
    [k1 k2] = size(abs1);
    
    winlen = 1000;   % window size
    fqlen = 1;       % freq. step
    entloopend = k2 / winlen;
    freqloopend = k1 / fqlen;
    
    Rhxabs = []; Rhyabs = []; Rhxyabs = [];     % preallocating output variables
    Rhycxabs = []; Rhxcyabs = [];
    Rixyabs = []; Ruxyabs = []; Ruyxabs = [];
    Rixynormabs = [];
    Rrelshanabs = [];
    
% ABS LOOP
    for entloop = 1:entloopend*2-1        % WINDOW LOOP
        inx1 = (entloop - 1) * winlen / 2 + 1;  % Note: overlaping windows!
        inx2 = inx1 + winlen - 1;
        
        RThxabs = []; RThyabs = []; RThxyabs = [];    % preallocating
        RThycxabs = []; RThxcyabs = [];
        RTixyabs = []; RTuxyabs = []; RTuyxabs = [];
        RTixynormabs = [];
        RTrelshanabs = [];
        
        for freqloop = 1:freqloopend       % FREQUENCY LOOP
            finx1 = (freqloop - 1) * fqlen + 1;
            finx2 = freqloop * fqlen;
            tfv1 = abs1(finx1:finx2,:);     % eeg wavelet "rows"
            tfv2 = abs2(finx1:finx2,:);     % unit wavelet "rows"
            y1 = abs1(finx1:finx2,inx1:inx2);   % eeg wavelet "cells"
            y2 = abs2(finx1:finx2,inx1:inx2);   % unit wavelet "cells"
            leny = length(y1);
            hb = fix(exp(0.626+0.4*log(leny-1)));   % number of bins
            [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entr2(tfv1,tfv2,y1,hb,y2,hb); % MAIN
            Ixynorm = Ixy / log(leny);
            relshan = Hy / log(hb);
            
            RThxabs = [RThxabs; Hx]; 
            RThyabs = [RThyabs; Hy]; 
            RThxyabs = [RThxyabs; Hxy]; 
            RThycxabs = [RThycxabs; Hycx]; 
            RThxcyabs = [RThxcyabs; Hxcy]; 
            RTixyabs = [RTixyabs; Ixy]; 
            RTuxyabs = [RTuxyabs; Uxy]; 
            RTuyxabs = [RTuyxabs; Uyx]; 
            RTixynormabs = [RTixynormabs; Ixynorm]; 
            Rrelshanabs = [Rrelshanabs; relshan];
        end         % end of freq. loop 
        
        Rhxabs = [Rhxabs RThxabs]; 
        Rhyabs = [Rhyabs RThyabs]; 
        Rhxyabs = [Rhxyabs RThxyabs]; 
        Rhycxabs = [Rhycxabs RThycxabs]; 
        Rhxcyabs = [Rhxcyabs RThxcyabs]; 
        Rixyabs = [Rixyabs RTixyabs]; 
        Ruxyabs = [Ruxyabs RTuxyabs]; 
        Ruyxabs = [Ruyxabs RTuyxabs]; 
        Rixynormabs = [Rixynormabs RTixynormabs]; 
        Rrelshanabs = [Rrelshanabs RTrelshanabs];
    end     % end of abs (window) loop
    
% Save
    fs = findstr(fnm,'_');
    fn = fnm(fs(2)+1:end-4);
    fnts = [wh_to fn '_RNDVSEEG'];
    save(fnts,'Rhxabs','Rhyabs','Rhxyabs','Rixyabs','Rixynormabs','Rrelshanabs','Rhxcyabs',...
        'Rhycxabs','Ruxyabs','Ruyxabs');
    clear abs1 abs2
    
    waitbar(i/sf)   % progress indicator
end     % end of cell cycle
close(wb)