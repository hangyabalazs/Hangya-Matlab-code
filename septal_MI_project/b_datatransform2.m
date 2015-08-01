function b_datatrasnsform2
%DATATRANSFORM2   Transforms sampling frequency to 10 kHz and converts data to usual form.
%   DATATRANSFORM2 resamples data registered on 16.667 kHz for the eeg and 8.333 kHz for
%   the unit on 10 kHz for both and saves two channel data (1.: eeg, 2.: unit).
%
%   See also DATATRANSFORM1, DATATRANSFORM3, DATATRANSFORM4 and RESAMPLE.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
global DATADIR
if isempty(DATADIR)
    clear glogal DATADIR;
    b_startup
    global DATADIR
end;
where1 = ['f:\raw_data\hcn\sampled_on_16_67_for_unit_and_8_33_for_eeg\temp\'];
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd(['f:\raw_data\hcn\']);  % here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    % progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    meret = size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end
    if size(data,2)==1
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
    
    datinx1 = 1;    % first point of the interval
    datinx2 = size(data,1);     % last point of the interval
    
    unit = data(:,2);
    unit = unit';
    eeg = data(:,1);
    eeg = eeg';
    
    eeg2 = eeg(1:round(2*length(eeg)/3));       % cut again
    unit2 = [eeg(round(2*length(eeg)/3)+1:end) unit];
    eeg = eeg2;
    unit = unit2;
    
    dt = 0.0001;
    time = [0:length(unit)-1] * dt; 
    xlimit = [min(time),max(time)];
    meret = size(data,1);
    mintafr = 1 / dt;
    
    global IN       % create global IN
    IN = cell(1,12);
    IN{1} = data;
    IN{2} = eeg;
    IN{3} = fname;
    IN{4} = where1;   % path name
    IN{5} = datinx1;
    IN{6} = datinx2;
    IN{7} = time;
    IN{8} = unit;
    IN{9} = dt;
    IN{10} = meret;
    IN{11} = mintafr;
    IN{12} = xlimit;
    
% Transform EEG
    eeg = resample(eeg,6,5);
    new_sampfreq = round(8333*6/5);
    
% Transform unit
%     [T,seglen] = b_thres3;      % thresholding
        
%     b_disc2(T,seglen);      % discrimination
    b_disc      % discrimination
    global DISC
    id = DISC{1};
    output = DISC{2};
    vdisc = DISC{3};
    
    dvr = 16667 / new_sampfreq;
    vdisc = vdisc / dvr;
    vdisc = round(vdisc);
    if vdisc(1) < 1     % corrigate rounding errors
        vdisc(1) = 1;
    end
    if vdisc(end) > length(eeg) - 1
        vdisc(end) = length(eeg) - 1;
    end
    unit = zeros(size(eeg));
    unit(vdisc) = 1;
    unit(vdisc+1) = -1;
    
% Transform data
    data = [eeg; unit]';
    
% Save
    eval(['save ' fname(1:6) ' data']);
    waitbar(o/sf)   % progress indicator
end 
close(wb);