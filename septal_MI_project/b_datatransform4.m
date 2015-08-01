function b_datatrasnsform4
%DATATRANSFORM4   Converts data to usual form.
%   DATATRANSFORM4 Saves two channel data (1.: eeg, 2.: unit).
%
%   See also DATATRANSFORM1, DATATRANSFORM2, DATATRANSFORM3 and RESAMPLE.

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
where1 = ['f:\raw_data\hcn\sampled_on_10_channel_no_3\'];
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
    
% Transform data
    data = data(:,1:2);
    
% Save
    eval(['save ' fname(1:6) ' data']);
    waitbar(o/sf)   % progress indicator
end 
close(wb);