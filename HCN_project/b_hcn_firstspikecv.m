function b_hcn_firstspikecv
%HCN_FIRSTSPIKECV   Compares first-spike variance and first-spike CV values.
%   HCN_FIRSTSPIKECV calculates and saves inter-first-spike interval variance
%   and coefficient of variation (CV) values for each analysed segment.
%
%   You are able to modify the input and the output directory through editing
%   the program code.
%
%   See also HCN_ANALYSIS and ITCLUST2.

% Input argument check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end
where = ['f:\raw_data\hcn\'];    % here are the data files
files = dir(where);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
    end
end
files2 = files2(2:end);
sf = length(files2);
mmm = pwd;
cd([DATAPATH,'HCN\Firstspikecv\']);  % here are the results

% Progress indicator
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Initialize output
mfspcv = [];
mfspvar = [];

% Import
for o = 1:sf    % "CELL CYCLE"
    fname = files2(o).name;
    ffnm = [where fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    data = load(ffnm);
    if isstruct(data)
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end
    if size(data,2)==1
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
    
% Computing the input variables
    meret = size(data,1);
    segmlen = 2000000;
    segno = floor(meret/segmlen) + 1;
    for seg = 1:segno   %"REGISTRATION SEGMENT CYCLE"
        datinx1 = (seg - 1) * segmlen + 1;      %first point of the interval
        datinx2 = min(seg*segmlen,meret);       %last point of the interval
        b_imp(fname,where,data,datinx1,datinx2);
        
% Thresholding
%         [T,seglen] = b_thres2;
%         cd thres
%         str = ['HCN_THRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'];
%         eval(['save ' str ' T seglen']);
%         cd ..

% Loading the saved threshold value
        directory = [DATAPATH,'HCN\Data\thresholds\'];
        str = [directory,'MEGATHRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'];
        load(str)
    
% Discrimination
%         b_disc2(T,seglen);
        b_disc(kuszob);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};
        
% Loading & analysing theta segments
        if segno == 1
            fln = ['THETA_SEGMENTS_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet\theta_segments',fln);
        else
            fln = ['THETA_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet\theta_segments_long',fln);
        end
        load(ff)
        segments = ThetaSegments;
        
        if ~isempty(segments)
            segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
        end                                         % shorter than 5 sec.
        
        th_index = size(segments,2);
        for t = 1:th_index      % "THETA SEGMENT CYCLE"
            
            segfirst = segments(1,t);
            seglast = segments(2,t);
            eeg_input = eeg(segfirst:seglast);
            unit_input = unit(segfirst:seglast);
            vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
            
            real_first_index = segfirst + datinx1 - 1;        % save segment type and bounderies
            real_last_index = seglast + datinx1 - 1;
            
            [mfspcv(end+1),mfspvar(end+1)] = hcn_main(eeg_input,unit_input,vdisc_input,...
                segfirst,seglast,filenam,real_first_index,real_last_index,'theta');
            
        end     % end of "theta segment cycle"
        
% Loading & analysing non-theta segments
        if segno == 1
            fln = ['NONTHETA_SEGMENTS_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet\nontheta_segments',fln);
        else
            fln = ['NONTHETA_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet\nontheta_segments_long',fln);
        end
        load(ff)
        segments = NonThetaSegments;
        
        if ~isempty(segments)
            segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
        end                                         % shorter than 5 sec.
        
        no_index = size(segments,2);
        for t = 1:no_index      % "NON-THETA SEGMENT CYCLE"
            
            segfirst = segments(1,t);
            seglast = segments(2,t);
            eeg_input = eeg(segfirst:seglast);
            unit_input = unit(segfirst:seglast);
            vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
            
            real_first_index = segfirst + datinx1 - 1;        % save segment type and bounderies
            real_last_index = seglast + datinx1 - 1;
            
            [mfspcv(end+1),mfspvar(end+1)] = hcn_main(eeg_input,unit_input,vdisc_input,...
                segfirst,seglast,filenam,real_first_index,real_last_index,'nontheta');
            
        end     % end of "non-theta segment cycle"
                
    end     % end of "registration segment cycle"
    
    waitbar(o/sf)   % progress indicator
end     % end of "cell cycle"
close(wb);   % close progress indicator

% Save output
mfspcv = mfspcv(find(mfspcv>=0));
mfspvar = mfspvar(find(mfspvar>=0));
save FirstSpikeVarCv mfspcv
save FirstSpikeVar mfspvar

% ----------------------------------------------------------------------------------
function segments = short_killer(segments)

% Skip short segments
int = segments;
int1 = int(1,:) + 200;    % leaving 200 - 200 ms edges of each segment
int2 = int(2,:) - 200;
difint = int2 - int1;
fd = find(difint<50000);         % leaving segments shorter than 5 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];

% ----------------------------------------------------------------------------------
function [mfspcv,mfspvar] = hcn_main(eeg,unit,vdisc,segfirst,seglast,filenam,real_first_index,real_last_index,st)

% Initialize output
mfspcv = -1;
mfspvar = -1;

% Ward's cluster analysis
enough = 1;
if length(vdisc) < 100
    enough = 0;
end

if enough
    [FirstSpikeCv,FirstSpikeVar,Burst] = itclust(vdisc,20,'ward');     %cluster analysis
    mfspcv = min(FirstSpikeCv(2:7));      % minimum of first-spike-cv (from 2 to 7 clusters)
    mfspvar = min(FirstSpikeVar(2:7));     % minimum of first-spike-var (from 2 to 7 clusters)
    
    Vdisc = vdisc;
    str = ['HCN_FIRST_SPIKE_CV_',filenam,'_',num2str(real_first_index),'_',num2str(real_last_index),'.mat'];
    eval(['save ' str ' FirstSpikeCv FirstSpikeVar Burst Vdisc']);

end

% ----------------------------------------------------------------------------------
function [firstspike_cv,firstspike_var,Burst] = itclust(vdisc,d,method)

% Input arguments check
error(nargchk(2,3,nargin));
if nargin == 2
    method = 'ward';    %set default method to 'ward';
end

% DIST, LINKAGE and COPHENET
ivs = diff(vdisc);
livs = length(ivs);
dmtx = zeros(livs,2);
iivs = ivs';
dmtx(:,1) = iivs;
dist = pdist2(dmtx);
links = linkage(dist,method);
Ccc = cophenet(links,dist);

% Allocating some matrices
firstspike_cv = zeros(1,d);
firstspike_var = zeros(1,d);

Burst = cell(1,d);  %allocating Burst cell

% CLUSTER
miniv = min(ivs);   %finding the cluster containing the smallest interval
for dec = 2:d
    c = cluster(links,dec);
    cell_clusters = cell(dec,1);
    doub_length = zeros(dec,1);
    for t = 1:dec
        cell_clusters{t} = find(c==t);
        doub_length(t) = length(cell_clusters{t});
        fnd = find(ivs(cell_clusters{t}) == miniv);
        if isempty(fnd) == 0,
            miniv_clus = t;
        end
    end
    
% Extraburst intervals
    ivss = ivs(cell_clusters{miniv_clus});  %creating extraburstiv
    extraburstiv = ivs(find(ivs>max(ivss)));
    
% Bursts
    liv = find(ivs<=max(ivss)); %finding the bursts
    liv1 = [-1 liv]; 
    liv = [liv 0];
    bliv = find(liv~=liv1+1);
    blivv = bliv(2:end)-1; 
    bliv(end) = [];
    burst = [liv(bliv);liv(blivv)+1];     %1st and last spikes of bursts
    fsp = vdisc(burst(1,:));              %localisation of 1st spikes of bursts
    lsp = vdisc(burst(2,:));              %localisation of last spikes of bursts
    diffburst = vdisc(burst(2,:)) - vdisc(burst(1,:));
    mdb = max(diffburst);
    
    Burst{dec} = burst; %filling Burst cell
    
% Inter-1st-sike intervals
    interfirstspike = diff(fsp);
    
% CV of the inter-1st-spike interval length
    if length(interfirstspike) ~= 0
        firstspike_cv(dec) = std(interfirstspike) / mean(interfirstspike);
    else firstspike_cv(dec) = NaN;
    end
    
% Variance of the normalized inter-1st-spike interval length
    if length(interfirstspike) ~= 0
        firstspike_var(dec) = var(interfirstspike / mean(interfirstspike));
    else firstspike_var(dec) = NaN;
    end
    
end