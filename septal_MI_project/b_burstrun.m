function [list,succr] = b_burstrun
%BURSTRUN    Saving firing parameters for burst selection.
%   [LIST, SR] = BURSTRUN calculates and saves parameters returned by 
%   ITCLUST, return plot data and K-means clustering of return plot data.
%   It works on theta and non-theta segments longer then 5 sec. (after
%   leaving 200-200 ms lag). It returns and saves success rates of K-means
%   clustering (SR) and list of analyzed files' names (LIST);
%
%   Note, that unit segment should contain at least 100 spikes for proper
%   clustering! Segments not satisfying this criterion are skipped.
%
%   User can change input and result directory, as well as directories
%   containing theta and non-theta segments through editing the program
%   code.
%
%   See also ITCLUST and KMEANS.

% Input argument check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'mi_zshift_data\discriminated_hc\'];    % here are the data files
resdir_theta = [DATAPATH,'Burst\Burst\theta_hc\'];  % here are the results
resdir_noth = [DATAPATH,'Burst\Burst\noth_hc\'];
othersdir = [DATAPATH,'Burst\Burst\others_hc\'];
thetadir = [DATAPATH,'Wavelet_hc\theta_segments\'];
nothdir = [DATAPATH,'Wavelet_hc\nontheta_segments\'];
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
    end
end
files2 = files2(2:end);
sf = length(files2);
list = cell(sf,1);      % initialize filename list
mmm = pwd;

% Progress indicator
wb = waitbar(0,'Running BURSTRUN...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
succr = cell(sf,2);    % initialize succes rate of K-means clustering
for o = 1:sf    % "CELL CYCLE"
    fname = files2(o).name;
    ffnm = [inpdir fname];
    filenam = [fname(1:3) '_' fname(5:6)];  % 'filenam' and 'filenam2' are different for 3ch data
    filenam2 = [fname(1:6)];
    list{o} = filenam2;
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    load(ffnm);
    
% Theta segments
    cd(resdir_theta)
    fln = ['THETA_SEGMENTS_',filenam];
    ff = fullfile(thetadir,fln);
    load(ff)
    segments = ThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

    th_index = size(segments,2);
    for t = 1:th_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
%         eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
        if length(vdisc_input) < 100
            continue
        end
        succr{o,1}(end+1) = burstrun_main(vdisc_input,filenam2,segfirst,seglast,'theta');
    end     % end of "theta segment cycle"
    
% Non-theta segments
    cd(resdir_noth)
    fln = ['NONTHETA_SEGMENTS_',filenam];
    ff = fullfile(nothdir,fln);
    load(ff)
    segments = NonThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end
    
    no_index = size(segments,2);
    for t = 1:no_index      % "NONTHETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
%         eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
        if length(vdisc_input) < 100
            continue
        end
        succr{o,2}(end+1) = burstrun_main(vdisc_input,filenam2,segfirst,seglast,'noth');
    end     % end of "non-theta segment cycle"
    waitbar(o/sf)
end     % end of "cell segment cycle"
close(wb)

% Save filename list and success rates
if ~b_isdir2('others2')
    mkdir others2
end
mm = pwd;
cd(othersdir)
save FileList list
save SuccessRate succr
cd(mm)



% ----------------------------------------------------------------------------------
function scr = burstrun_main(vdisc,filenam,segfirst,seglast,str)

% Iterative Cluster Analysis
d = 20;     % maximal number of clusters
method = 'Ward';
[IntraBurstIvCv,ExtraBurstIvCv,InterBurstIvCv,FirstSpikeCv,AllFirstSpikeCv,BurstLengthCv,...
    Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,...
    Burst,Ccc] = b_itclust2(vdisc,d,method);     %cluster analysis and burst parameters

% Return plot
ivs = diff(vdisc);
livs = length(ivs);
ReturnPlotXData = ivs(1:livs-1);
ReturnPlotYData = ivs(2:livs);

% K-means clustering on the returnplot data
[IdX,Centroid,SumDist,dd,totsumd] = b_kmeans2([ReturnPlotXData; ReturnPlotYData]',3,'replicates',10);
scr = length(find(totsumd==min(totsumd))) / length(totsumd);

% Vdisc
Vdisc = vdisc;

% Save
fn = [filenam '_' str '_' num2str(segfirst) '_' num2str(seglast) '_BURST'];
str1 = ['save ' fn ' IntraBurstIvCv ExtraBurstIvCv InterBurstIvCv'];
str2 = [' FirstSpikeCv AllFirstSpikeCv BurstLengthCv Burstiness'];
str3 = [' IntraBurstFrequency IntraBurstSpikeNumber BurstLength'];
str4 = [' BurstFrequency Burst Ccc ReturnPlotXData ReturnPlotYData'];
str5 = [' IdX Centroid SumDist'];
str6 = [' Vdisc'];
str = [str1 str2 str3 str4 str5 str6];
eval(str)



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