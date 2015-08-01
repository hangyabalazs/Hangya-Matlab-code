function b_thresrun2
%THRESRUN2   Runs automatic thresholder on a sequence of files.
%   THRESRUN2 uses two directories: one for the input files and another for the results;
%   you are able to modify these directories through editing the program code.
%
%   THRESRUN2 saves threshold values in the results' directory. After running THRESRUN2
%   you are able to create new type data files (containing eeg and vdisc) using THRES_GUI.
%
%   See also THRES5, DISC2, DISCRUN and THRES_GUI.

% Input arguments check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files\temp\'];    %Here are the data files
resdir = [DATAPATH,'Raphe\raphe_juxta\Thres\'];            %Here are the results
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;
cd(resdir);

% Import
wb = waitbar(0,'Running THRESRUN2...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    fname = files(o).name;
    filenam = fname(1:6);
    ffnm = [inpdir fname];
    try
        load(ffnm)
        unit = Unit.values';
        eeg = hEEG.values';
    catch
        data = b_load_data(ffnm);
        unit = data(:,2)';
        eeg = data(:,1)';
    end
    dt = 0.0001;
    time = [0:length(unit)-1] * dt;
        
% Thresholding
    close all
    [T,seglen] = thres(unit,time);
    
% Discrimination
    vdisc = disc(unit,T,seglen);
    str = ['THRES_',fname];
    save(str,'T','seglen','vdisc');
    
% Saving
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);

% -------------------------------------------------------------------------
function vdisc = disc(unit,kuszob2,seglen)

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
    [maxe,maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end

% -------------------------------------------------------------------------
function [T,seglen] = thres(unit,time)

% Assign segment length
seglen = 100000;

% Threshold pseudounit
if length(find(unit==max(unit))) > 15 & isequal(max(unit),1)
    T = [];
    ind2 = 0;
    lenu = length(unit);
    while ind2 < lenu
        ind2 = ind2 + seglen;
        T(end+1) = 0.25;
    end
    return
end

% Checking if base line is stable
d = 10000;
dd = 1;
mn = [];
while dd < time(end)
    mn(end+1) = mean(unit(dd:dd+d));
    dd = dd + d;
end
mn(end+1) = mean(unit(dd:end));
if std(mn) > 0.1
    str = [fname(1:6) ' : Instable base line.'];
    disp(str)
end

% Plot unit
% H = figure;
% plot(time,unit)

% Segmenting the unit
MS = [];
MI = [];
T = [];
ind2 = 0;
lenu = length(unit);
while ind2 < lenu
    if lenu - ind2 < 30000      % too short segment can cause false thresold
        v = T(end);
        T(end+1) = v;
        break
    end
    ind1 = ind2 + 1;
    ind2 = ind2 + seglen;
    unitsegment = unit(ind1:min(ind2,lenu));
    
% Action potential amplitudes
    unit2 = unitsegment;
    for k = 1:20
        m(k) = max(unit2);
        fnd = find(unit2==m(k));
        unit2(fnd(1)) = [];
    end
    MS(end+1) = min(m(16:20));
    
% Noise
    z = mean(unitsegment);
    m = max(unitsegment);
    npc = 1;
    nsup = 0.98;
    ninf = 0.97;
    next = 0;
    while (npc > nsup | npc < ninf) & next <= 200
        if npc > ninf
            m = m - (m - z) / 2;
        else
            m = m + (m - z) / 2;
        end
        f = find(unitsegment<m);
        npc = length(f) / length(unitsegment);
        next = next + 1;
        if next == 201
            str = [fname(1:6) ' : Maximum number of iterations reached in noise limit search.'];
            disp(str)
        end
    end
    MI(end+1) = m;
    
% Threshold
    T(end+1) = MI(end) + (MS(end) - MI(end)) / 4;

% Plot threshold
%     line([time(ind1) time(min(ind2,lenu))],[MI(end) MI(end)],'Color','r')
%     line([time(ind1) time(min(ind2,lenu))],[MS(end) MS(end)],'Color','r')
%     line([time(ind1) time(min(ind2,lenu))],[T(end) T(end)],'Color','g')
end
% t = ['THRESHOLD ',fname(1:3),' ',fname(5:6),' ',num2str(datinx1),' ',num2str(datinx2)];
% title(t);