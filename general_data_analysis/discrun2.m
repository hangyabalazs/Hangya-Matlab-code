function discrun2
%DISCRUN2   Runs automatic thresholder and discriminator on a sequence of files.
%   DISCRUN2 uses two directories: one for the input files and another for
%   the results; you are able to modify these directories through editing
%   the program code.
%
%   DISCRUN2 saves threshold values and 'vdisc' in the results' directory.
%   After running DISCRUN2 you are able to create new data files
%   (containing eeg and vdisc) using THRES_GUI.
%
%   See also THRES4, DISC2, THRESRUN2 and THRES_GUI.

% Input arguments check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
global DATADIR
% inpdir = [DATADIR 'PViktor\raw\'];
inpdir = [DATADIR 'PViktor\non-glycinergic\'];
resdir = [DATAPATH 'PViktor\thres\'];            %Here are the results
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;
cd(resdir);

% Import
wb = waitbar(0,'Running DISCRUN...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    fname = files(o).name;
    ffnm = [inpdir fname];
    data = load(ffnm);
    unit = data.unit.values';
    sr = round(1/(data.unit.times(2)-data.unit.times(1)));
    
    % Thresholding
    [T,seglen] = thres5(unit,sr);
    
    % Discrimination
    vdisc = ldisc(unit,T,seglen); %#ok<NASGU>
    
    % Saving
    str = ['THRES_',fname(1:end-4),'.mat'];
    eval(['save ' str ' T seglen vdisc']);
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
cd(mm);

% -------------------------------------------------------------------------
function vdisc = ldisc(unit,kuszob2,seglen)
%DISC   Discriminator for raw unit data.

% Threshold
ind2 = 0;
lenu = length(unit);
next = 1;
kuszob = [];
while ind2 < lenu
    ind1 = ind2 + 1;
    ind2 = ind2 + seglen;
    kuszob(ind1:min(ind2,lenu)) = kuszob2(next);
    next = next + 1;
end

% Discriminating
disc = find(unit>=kuszob); 
discl = length(disc);
disc = [disc; unit(disc(1:discl))];
disc(1,discl+1) = length(unit) + 2;
dif = diff(disc(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j = 1:length(difn1) - 1
    [~, maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end