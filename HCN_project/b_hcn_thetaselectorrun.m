function b_hcn_thetaselectorrun
%HCN_THETASELECTORRUN   Runs THETASELECTOR3 and THETASELECTOR_BETA3 on a sequence of files.
%   HCN_THETASELECTORRUN uses two directories: one for the data files and one for the 
%   results. You are able to modify these directories through editing the program code.
%
%   The function saves the output matrix (OM) in the 'matrix' subdirectory. OM is a 
%   3-N matrix, which contains THETASELECTOR3 output matrix (maximum locations and maximum
%   values) in the first and second row, and THETASELECTOR_BETA3 output matrix (theta-delta
%   ratio) in the third row. THETASELECTOR3 output figure (wavelet with maximum locations) 
%   gets saved in 'jpg' subdirectory. Subdirectories should preexist in the results'
%   directory.
%
%   You have to use HCN_THETASELECTORRUN_LONG for registrations longer than 200 seconds.
%   HCN_THETASELECTORRUN saves a binary file called 'long.mat', which contains the 
%   identifiers of the long registrations.
%
%   See also THETASELECTOR3, THETASELECTOR_BETA3, THETASELECTORRUN, THETASELECTORRUN_LONG
%   and HCN_THETASELECTORRUN_LONG.

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
% where1 = DATADIR;    %Here are the data files
where1 = ['f:\raw_data\hcn\all\'];
files = dir(where1);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
    end
end
files2 = files2(2:end);
sf = length(files2);
mm = pwd;
cd([DATAPATH, 'HCN\Wavelet2\thetaselection\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf,
    fname = files2(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    meret = size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2)==1,
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
    
% Computing the input variables
    datinx1 = 1; %first point of the interval
    datinx2 = size(data,1); %last point of the interval
    if datinx2 > 2000000
        global LONG
        LONG{end+1} = fname(1:6);
        waitbar(o/sf)   %Progress indicator
        continue
    end
    
    dt = 0.0001;
    mintafr = 1 / dt;
    unit = data(datinx1:datinx2,2);
    unit = unit';
    eeg = data(datinx1:datinx2,1);
    eeg = eeg';
    time = [0:length(unit)-1] * dt; 
    xlimit = [min(time),max(time)];
    global IN
    IN = cell(1,12);
    IN{1} = data;
    IN{2} = eeg;
    IN{3} = fname;
    IN{4} = where1;     % pathname
    IN{5} = datinx1;
    IN{6} = datinx2;
    IN{7} = time;
    IN{8} = unit;
    IN{9} = dt;
    IN{10} = meret;
    IN{11} = mintafr;
    IN{12} = xlimit;
    
% Free memory
    save temp o sf wb files2 where1
    clear
    clear function
    
% Theta selection
    [wave,f,newstep] = b_waveletcall_for_applications;
    save newstep newstep
    save scalefreq f
    
    sw2 = size(wave,2);
    pieceno = 5;
    segm = fix(sw2/pieceno);
    power = abs(wave) .^ 2;
%     power = [];
%     while ~isempty(wave)
%         index1 = 1;
%         index2 = min(segm,size(wave,2));
%         wavefrag = wave(:,index1:index2);
%         powerfrag = (abs(wavefrag)) .^ 2;
%         clear wavefrag
%         wave(:,index1:index2) = [];
%         power = [power powerfrag];
%     end
    clear wave
    
    [H,OM] = b_thetaselector3(power,f,newstep);
    ratio = b_thetaselector_beta3(power,f,newstep);
    
% Saving
    Out = zeros(3,sw2);
    Out(1:2,:) = OM;
    Out(3,:) = ratio;
    
    global IN
    fname = IN{3};
    pont = findstr(fname,'.');
    filenam = fname(1:pont(1)-1);
    fig = gcf;
    ax = findobj(fig,'Type','axes');
    axes(ax(2))
    title(['EEG ',filenam(1:3),' ',filenam(5:6)]);
    eval(['saveas(H,''THETA_SELECT_',filenam(1:6),'.jpg'')']);
    eval(['save(''THETA_SELECT_',filenam(1:6),'.mat'',''Out'')']);
    load temp
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
warning off
delete temp.mat
warning on
close all

global LONG
long = LONG;
save long long
clear global LONG