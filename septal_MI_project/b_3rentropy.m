function b_3rentropy
%3RENTROPY   Entropy calculation.
%   3RENTROPY calculates entropies for 5 sec. long, 80 % overlaping windows. The concidered
%   distributions are interspike interval distribution and eeg amplitude distribution. It 
%   calls ENTROPY for the raw entropy calculation. Outputs:
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
%   See also ENTROPY and RENTROPY.

% Input argument check
error(nargchk(0,0,nargin))

% Define directories
global DATADIR
global DATAPATH
wh_from = [DATAPATH 'Data\Analysenow3\'];
wh_to = [DATAPATH 'Entropy3\raw\real\'];
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
wb = waitbar(0,'Running 3RENTROPY...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
for i = 1:sf    % CELL CYCLE
    fnm = files2(i).name;
    ffnm = [wh_from fnm];
    load(ffnm);
    
    abs1 = eeg;
    isi = diff(vdisc);
    abs2 = zeros(1,length(eeg));
    for i = 1:length(vdisc)-1
        abs2(vdisc(i):vdisc(i+1)) = isi(i);
    end
    [k1 k2] = size(abs1);
    
    winlen = 50000;   % window size
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
        
        y1 = abs1(inx1:inx2);   % eeg
        y2 = abs2(inx1:inx2);   % unit
        leny = numel(y1);
        hb = fix(exp(0.626+0.4*log(leny-1)));   % number of bins
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entr2(abs1,abs2,y1,hb,y2,hb); % MAIN
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
    fn = fnm(1:6);
    fnts = [wh_to fn '_RENTROPY'];
    save(fnts,'Rhxabs','Rhyabs','Rhxyabs','Rixyabs','Rixynormabs','Rrelshanyabs','Rrelshanxabs',...
        'Rhxcyabs','Rhycxabs','Ruxyabs','Ruyxabs');
    clear abs1 abs2
    
    waitbar(i/sf)   % progress indicator
end     % end of cell cycle
close(wb)