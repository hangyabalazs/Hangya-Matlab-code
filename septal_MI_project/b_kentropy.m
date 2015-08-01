function b_kentropy
%KENTROPY   Entropy calculation.
%   KENTROPY calculates entropies for whole data segments. The concidered distributions are
%   interspike interval distribution and eeg amplitude distribution. It calls ENTROPY for the
%   raw entropy calculation. Outputs:
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
%   See also ENTROPY.

% Input argument check
error(nargchk(0,0,nargin))

% Define directories
global DATADIR
global DATAPATH
wh_from = [DATAPATH 'Data\Analysenow3\'];
wh_to = [DATAPATH 'Entropy2\whole\real\'];
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
wb = waitbar(0,'Running KENTROPY...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
for i = 1:sf    % CELL CYCLE
    fnm = files2(i).name;
    ffnm = [wh_from fnm];
    load(ffnm);
    
    y1 = eeg;
    isi = diff(vdisc);
    y2 = zeros(1,length(eeg));
    for i = 1:length(vdisc)-1
        y2(vdisc(i):vdisc(i+1)) = isi(i);
    end
    leny = numel(y1);
    hb = fix(exp(0.626+0.4*log(leny-1)));   % number of bins
    
    [hx,hy,jh,Rhxabs,Rhyabs,Rhxyabs,Rhxcyabs,Rhycxabs,Rixyabs,Ruxyabs,Ruyxabs] = b_entropy(y1,y2);
    Rixynormabs = Rixyabs / log(leny);
    Rrelshanyabs = Rhyabs / log(hb);
    Rrelshanxabs = Rhxabs / log(hb);
    
% Save
    fn = fnm(1:6);
    fnts = [wh_to fn '_WENTROPY'];
    save(fnts,'Rhxabs','Rhyabs','Rhxyabs','Rixyabs','Rixynormabs','Rrelshanyabs','Rrelshanxabs',...
        'Rhxcyabs','Rhycxabs','Ruxyabs','Ruyxabs');
    clear abs1 abs2
    
    waitbar(i/sf)   % progress indicator
end     % end of cell cycle
close(wb)