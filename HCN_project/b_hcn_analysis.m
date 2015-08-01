function b_hcn_analysis
%HCN_ANALYSIS   Analysis of HCN cell data.
%   This function uses two main directories - one for the data files and one for the results. 
%   You have to specify these directories in the program code. (Subdirectories in the
%   results' directory are created automatically.) Some other directories are also used:
%   directories containing the segment boundaries and the thresholds.
%
%   HCN_ANALYSIS analyses theta, nontheta and sharp wave segments seperately. The segment
%   boundaries are loaded from the appropriate directories. (For details on segment selection,
%   see HCN_THETASELECTORRUN, THETASELECTOR3, THETA and NONTHETA.) The segments have to be 
%   longer than 5 sec. and contain more than 100 spikes.
%
%   The thresholds for unit discrimination are former saved and here loaded: see HCN_MEGATHRES
%   for details. Discrimination is done by DISC.
%
%   The unit analysis is based on Ward's cluster analysis. HCN_ANALYSIS creates different
%   number of clusters from 2 to 20 and computes intraburst interval, interburst interval,
%   "inter-fist-spike interval", "all-inter-first-spike-interval" (distances of first spikes 
%   including single spikes as well), extraburst interval length CV and burst length CV 
%   (coefficient of variation: std(X) / mean(X)) at all number of clusters.
%   Burst parameters are also computed: burstiness (number of intraburst spikes / number of
%   all spikes), average intraburst frequency (calculated for each burst and averaged),
%   average intraburst spike number, average burst length, frequency of burst first spikes
%   (called 'burst frequency'). If the segment is nonbursty, all these parameters are assigned
%   as NaN. See also ITCLUST.
%   HCN_ANALYSIS applies the following burst definition: intraburst intervals are the interspike
%   intervals containted by the interspike interval cluster with the shortest interval.
%   The appropriate cluster number is selected and saved via HCN_GUI (semiautomatic burst
%   selection).
%
%   Average frequency is calculated for all type of segments.
%
%   Unit phase values are calculated for theta segments from eeg-unit crosswavelet phase.
%   Angle phase histogram and circular statistics (circular mean, circular standard deviation,
%   mean vector length) are calculated and saved. For a detailed description of wavelet analysis,
%   see WAVELET_NEW3.
%
%   Relative theta power of unit is calculated from unit wavelet power.
%
%   Analysis:
%       I. Segment selection
%               A) theta
%               B) nontheta
%               C) sharp waves
%       II. Firing pattern analysis
%               A) firing rate
%               B) relative theta power
%               C) burst analysis
%                   0. burst selection
%                   1. burstiness
%                   2. intraburst frequency
%                   3. intraburst spike number
%                   4. burst frequency
%                   5. burst length
%               D) variability analysis
%                   1. intraburst frequency CV
%                   2. interburst frequency CV
%                   3. extraburst frequency CV
%                   4. firstspike frequency CV
%                   5. allfirstspike frequency CV
%                   6. burst length CV
%       III. Phase analysis
%               A) phase angles
%               B) phase histogram
%               C) circular statistics
%                   1. circular mean
%                   2. circular standard deviation
%                   3. mean vector length
%
%   Output:
%       1. hcn_analysis.txt: contain the firing pattern analysis data
%       2. "HCN_ICA" mat files: contain the iterative cluster analysis result and vdisc
%               (output of unit discrimination)
%       3. "HCN_BURST" mat files: contain the burst boundaries and vdisc
%       4. "PHASE_DIFF" mat files: contain the phase angle matrix
%       5. "PHASE_HIST" mat files: contain phase angle histogram, circular statistics
%               data and vdisc
%
%   See also HCN_MEGATHRES, HCN_BURST, HCN_THETASELECTORRUN, THETASELECTOR3 and ITCLUST.


% Input argument check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end
where = ['f:\raw_data\hcn\all\'];    % here are the data files
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

cd([DATAPATH,'HCN\Analysis5\']);  % here are the results
if ~b_isdir2('thres')
    mkdir thres
end
if ~b_isdir2('wavelet_settings')
    mkdir wavelet_settings
end

npt = input('Discard existing content or append data while writing hcn_analysis.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('hcn_analysis.txt','w');
elseif npt == 'a',
    fid = fopen('hcn_analysis.txt','a');
else
    error('Unexpected answer for input.')
end;

% Progress indicator
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
list = {};      % initialize list of analysed files
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
    if size(data,2) == 1
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
            ff = fullfile(DATAPATH,'HCN\Wavelet2\theta_segments',fln);
        else
            fln = ['THETA_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet2\theta_segments_long',fln);
        end
        load(ff)
        segments = ThetaSegments;
        
        if ~isempty(segments)
            segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
        end
        
        th_index = size(segments,2);
        for t = 1:th_index      % "THETA SEGMENT CYCLE"
            
            segfirst = segments(1,t);
            seglast = segments(2,t);
            eeg_input = eeg(segfirst:seglast);
            unit_input = unit(segfirst:seglast);
            vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
            
            real_first_index = segfirst + datinx1 - 1;        % save segment type and bounderies
            real_last_index = seglast + datinx1 - 1;
            fprintf(fid,'%s %s %f %f',filenam,'theta',real_first_index,real_last_index);
            
            hcn_main1(vdisc_input,segfirst,seglast,fid)
            if length(vdisc_input) >= 2
                hcn_main2(eeg_input,unit_input,vdisc_input,segfirst,seglast,fid,filenam,real_first_index,real_last_index,'theta')
            end
            
            fprintf(fid,'\n');
            
        end     % end of "theta segment cycle"
        
% Loading & analysing non-theta segments
        if segno == 1
            fln = ['NONTHETA_SEGMENTS_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet2\nontheta_segments',fln);
        else
            fln = ['NONTHETA_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet2\nontheta_segments_long',fln);
        end
        load(ff)
        segments = NonThetaSegments;
        
        if ~isempty(segments)
            segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
        end
        
        no_index = size(segments,2);
        for t = 1:no_index      % "NON-THETA SEGMENT CYCLE"
            
            segfirst = segments(1,t);
            seglast = segments(2,t);
            eeg_input = eeg(segfirst:seglast);
            unit_input = unit(segfirst:seglast);
            vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
            
            real_first_index = segfirst + datinx1 - 1;        % save segment type and bounderies
            real_last_index = seglast + datinx1 - 1;
            fprintf(fid,'%s %s %f %f',filenam,'nontheta',real_first_index,real_last_index);
            
            hcn_main1(vdisc_input,segfirst,seglast,fid)
            if length(vdisc_input) >= 2
                hcn_main2(eeg_input,unit_input,vdisc_input,segfirst,seglast,fid,filenam,real_first_index,real_last_index,'nontheta')
            end
            
            fprintf(fid,'\n');
            
        end     % end of "non-theta segment cycle"
        
% Loading & analysing sharp wave segments
        if segno == 1
            fln = ['SHARPWAVE_SEGMENTS_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet2\sharpwave_segments',fln);
        else
            fln = ['SHARPWAVE_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'HCN\Wavelet2\sharpwave_segments_long',fln);
        end
        load(ff)
        segments = SharpWaveSegments;
        
        sp_index = size(segments,2);
        for t = 1:sp_index      % "SHARP WAVE SEGMENT CYCLE"
            
            segfirst = segments(1,t);
            seglast = segments(2,t);
            eeg_input = eeg(segfirst:seglast);
            unit_input = unit(segfirst:seglast);
            
            real_first_index = segfirst + datinx1 - 1;        % save segment type and bounderies
            real_last_index = seglast + datinx1 - 1;
            fprintf(fid,'%s %s %f %f',filenam,'sharpwave',real_first_index,real_last_index);
            
            vdisc_input = vdisc(find(vdisc>segfirst&vdisc<seglast));
            hcn_main1(vdisc_input,segfirst,seglast,fid)
            
            fprintf(fid,'\n');
            
        end     % end of "sharp wave segment cycle"
        
    end     % end of "registration segment cycle"
    
    list{end+1} = [filenam];
    waitbar(o/sf)   % progress indicator
end     % end of "cell cycle"
close(wb);   % close progress indicator
fclose(fid);    % close text file id

% Save list of analysed files
list = list';
save list list
cd(mmm);

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
function hcn_main1(vdisc,segfirst,seglast,fid)

% Average frequency
AverageFrequency = length(vdisc) * 10000 / (seglast - segfirst);    % average frequency in Herz

% Save
fprintf(fid,'%s %f',' average_frequency',AverageFrequency);

% ----------------------------------------------------------------------------------
function hcn_main2(eeg,unit,vdisc,segfirst,seglast,fid,filenam,real_first_index,real_last_index,st)

% Ward's cluster analysis
enough = 1;
if length(vdisc) < 100
    enough = 0;
end

if enough
    global DATAPATH
    bwhere = [DATAPATH 'HCN\Burst\bursty\'];
%     nwhere = [DATAPATH 'HCN\Busrst\nonbursty\'];
    fby = dir(bwhere);
%     fnb = dir(nwhere);
    ns_long = cell(1,length(fby));
    ns_short = cell(1,length(fby));
    for i = 3:length(fby)
        ns_long{i} = fby(i).name;
        fs = findstr(ns_long{i},'_');
        ns_short{i} = ns_long{i}(fs(2)+1:end);
        ns_shorter{i} = ns_short{i}(1:end-4);
    end
    fname = [filenam '_' num2str(segfirst) '_' num2str(seglast)];
    scp = strcmp(fname,ns_shorter);
    fnd = find(scp);
    if isempty(fnd)
        isbursty = 0;
    else
        ff = [bwhere 'OPTIMAL_CLUSNO_' fname];
        load(ff)
        isbursty = 1;
    end
end

if enough & isbursty
    [IntraBurstIvCv,ExtraBurstIvCv,InterBurstIvCv,FirstSpikeCv,AllFirstSpikeCv,BurstLengthCv...
            Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,...
            Burst,Ccc] = b_itclust2(vdisc,20,'ward');     %cluster analysis and burst parameters
    Vdisc = vdisc;
    str = ['HCN_ICA_',filenam,'_',num2str(real_first_index),'_',num2str(real_last_index),'.mat'];
    eval(['save ' str ' IntraBurstIvCv ExtraBurstIvCv InterBurstIvCv FirstSpikeCv AllFirstSpikeCv BurstLengthCv Burst Vdisc']);
    
% Burst analysis
    clusinx = opt_clusno;
    if isbursty     % bursts: cluster analysis result at cluter no. of minimal first-spike-cv
        IntraBurstIvCv = IntraBurstIvCv(clusinx);   % variance analysis
        ExtraBurstIvCv = ExtraBurstIvCv(clusinx);
        InterBurstIvCv = InterBurstIvCv(clusinx);
        FirstSpikeCv = FirstSpikeCv(clusinx);
        AllFirstSpikeCv = AllFirstSpikeCv(clusinx);
        BurstLengthCv = BurstLengthCv(clusinx);
        
        Burstiness = Burstiness(clusinx);           % burst parameters
        IntraBurstFrequency = IntraBurstFrequency(clusinx);
        IntraBurstSpikeNumber = IntraBurstSpikeNumber(clusinx);
        BurstLength = BurstLength(clusinx);
        BurstFrequency = BurstFrequency(clusinx);
        
    else
        IntraBurstIvCv = NaN;   % variance analysis
        ExtraBurstIvCv = NaN;
        InterBurstIvCv = NaN;
        FirstSpikeCv = NaN;
        AllFirstSpikeCv = NaN;
        BurstLengthCv = NaN;
        
        Burstiness = NaN;           % burst parameters
        IntraBurstFrequency = NaN;
        IntraBurstSpikeNumber = NaN;
        BurstLength = NaN;
        BurstFrequency = NaN;
    end
    
    Burst = Burst{clusinx};                     % burst bounaries
    Vdisc = vdisc;
    str = ['save(''HCN_BURST_',filenam,'_',num2str(real_first_index),'_',num2str(real_last_index),'.mat'',''Burst'',''Vdisc'')'];
    eval(str);     % save 
        
    fprintf(fid,'%s %f',' intraburst_iv_cv',IntraBurstIvCv);     % save
    fprintf(fid,'%s %f',' extraburst_iv_cv',ExtraBurstIvCv);
    fprintf(fid,'%s %f',' interburst_iv_cv',InterBurstIvCv);
    fprintf(fid,'%s %f',' first_spike_cv',FirstSpikeCv);
    fprintf(fid,'%s %f',' all_first_spike_cv',AllFirstSpikeCv);
    fprintf(fid,'%s %f',' burst_length_cv',BurstLengthCv);
    fprintf(fid,'%s %f',' burstiness',Burstiness);
    fprintf(fid,'%s %f',' intraburst_frequency',IntraBurstFrequency);
    fprintf(fid,'%s %f',' intraburst_spike_number',IntraBurstSpikeNumber);
    fprintf(fid,'%s %f',' burst_length',BurstLength);
    fprintf(fid,'%s %f',' burst_frequency',BurstFrequency);
end

% Wavelet
[wave_eeg,wave_unit,f,newstep] = waveletcall(eeg,unit,vdisc-segfirst);      % wavelet transformation

ScaleVector = f(1:size(wave_eeg,1));     % save settings
Newstep = newstep;
cd wavelet_settings
save newstep Newstep
save f ScaleVector
cd ..

% Relative theta power
fnd = find(f>6);    % theta band bounderies
pwind1 = fnd(end);
fnd = find(f<3);
pwind2 = fnd(1);

unit_power = (abs(wave_unit)) .^ 2;     % unit power
thetapower = unit_power(pwind1:pwind2,:);    % theta power
sumthetapower = sum(thetapower);
clear thetapower
allpower = sum(unit_power);     % all power
clear unit_power
rtp = sumthetapower ./ allpower;
RelativeThetaPower = mean(rtp);

fprintf(fid,'%s %f',' relative_theta_power',RelativeThetaPower);      % save

% Crosswavelet
if strcmp(st,'theta')
    crosswave = wave_eeg .* conj(wave_unit);    % calculate crosswavelet
    clear wave_eeg wave_unit
    
    sw1 = size(crosswave,1);    % calculate crosswavelet power and phase
    sw2 = size(crosswave,2);
    pieceno = 5;
    segm = fix(sw2/pieceno);
    crosspower = [];
    crossphase = [];
    while ~isempty(crosswave)
        index1 = 1;
        index2 = min(segm,size(crosswave,2));
        wavefrag = crosswave(:,index1:index2);
        powerfrag = (abs(wavefrag)) .^ 2;
        phasefrag = angle(wavefrag);
        clear wavefrag
        crosswave(:,index1:index2) = [];
        crosspower = [crosspower powerfrag];
        crossphase = [crossphase phasefrag];
    end
    clear crosswave
    
% Phase angles
    maxes = max(crosspower);    % find maximums in crosspower
    maxloc = zeros(1,sw2);
    angs = [];
    for i = 1:sw2
        maxloc(i) = find(crosspower(:,i)==maxes(i));
        if maxloc(i) > pwind1 & maxloc(i) < pwind2
            angs(end+1) = crossphase(maxloc(i),i);
        end
    end
    ThetaPhaseDifference = angs / pi * 180;
    
    edges = [-pi : pi/30 : pi];
    AngleHistogram = histc(angs,edges);
    
    CircularMean = b_circular_mean2(angs);
    CircularStd = b_circular_std(angs);
    MeanVectorLength = b_mvl(angs);
    
    str = ['save(''HCN_PHASEDIFF_',filenam,'_',num2str(real_first_index),'_',num2str(real_last_index),...
            '.mat'',''ThetaPhaseDifference'')'];       % save
    eval(str);
    str = ['save(''HCN_PHASEHIST_',filenam,'_',num2str(real_first_index),'_',num2str(real_last_index),...
            '.mat'',''AngleHistogram'',''CircularMean'',''CircularStd'',''MeanVectorLength'',''vdisc'')'];
    eval(str);
end

% ----------------------------------------------------------------------------------
function [wave_eeg,wave_unit,f,newstep] = waveletcall(eeg,unit,vdisc)

% Create variables for wavelet
newstep = 25;
resamp = 10000/newstep;

sst_eeg = eeg(1:newstep:end);   % downsample eeg
sst_unit = b_sincconv2(unit,vdisc,newstep);    % unit sinc convolution

% Standardization
sst_eeg = (sst_eeg - mean(sst_eeg)) / std(sst_eeg);
sst_unit = (sst_unit - mean(sst_unit)) / std(sst_unit);
if abs(length(sst_unit)-length(sst_eeg)) > 5
    error('Different eeg and unit wavelet input length.')
end
minlen = min(length(sst_eeg),length(sst_unit));     % length sinchronization
if minlen == length(sst_eeg)
    sst_unit = sst_unit(:,1:minlen);
else
    sst_eeg = sst_eeg(:,1:minlen);
end

% Wavelet transformation
dt = 1 / resamp;    %resample on 400 Hz
time = [0:length(sst_unit)-1] * dt + 0;
n = length(sst_unit); 
pad = 1;      
dj = 0.04;
s0 = 2 * dt;
j1 = ((1 / dj) * log2(n/2)) * 2;
j1 = ceil(j1);
j = (0:j1);
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod; 
lag1 = 0.72;  
mother = 'Morlet';
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
[wave_eeg,period,scale,coi] = b_wavelet_new3(sst_eeg,dt,pad,dj,s0,j1,mother,param,mis);
[wave_unit,period,scale,coi] = b_wavelet_new3(sst_unit,dt,pad,dj,s0,j1,mother,param,mis);