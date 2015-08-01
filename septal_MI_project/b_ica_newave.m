function b_ica_newave(d)
%ICA_NEWAVE  Iterative Cluster Analysis and wavelet average vector calculation.
%   This function uses two directories - one for the data files and one for the results. 
%   You have to specify these directories in the program code.
%
%   ICA_NEWAVE(D) creates different number of clusters from 2 to D and computes intraburst 
%   interval, interburst interval, "inter-fist-spike interval", "all-inter-first-spike-
%   interval" (distances of first spikes including single spikes as well) and extraburst
%   interval length CV (coefficient of variation: std(X) / mean(X)) at all number of clusters.
%   ICA_NEWAVE, by itself, uses D = 20 default value.
%   Burst parameters are also computed: burstiness (number of intraburst spikes / number of
%   all spikes), average intraburst frequency (calculated for each burst and averaged),
%   average intraburst spike number, average burst length, frequency of burst first spikes
%   (called 'burst frequency').
%   Clustering methods 'Ward', 'single' and 'average' were used (see LINKAGE for details).
%
%   ICA_NEWAVE, unlike ICA_BETA creates only matricies, not figures. You can display the figures
%   using ICA_GUI3.
%
%   ICA_NEWAVE applies the following burst definition: intraburst intervals are the interspike
%   intervals containted by the interspike interval cluster with the shortest interval.
%
%   In contrast to ICA_WAVE, what works on one predefinied interval and loads previously
%   saved threshold files, ICA_NEWAWE works on the whole registration (spliting it at 200 sec.
%   if too long) and uses automatic thresholding (see THRES for details).
%
%   The program runs on theta segments -  selection is based on THETASELECTOR. The segments
%   have to be longer than 5 sec. and contain more than 100 spikes.
%
%   Wavelet spectrum is calculated using WAVELET function.
%
%   Outputs: 
%       1. Wavelet average magnitude vectors ('Wavevec') contain the mean along scale of 
%          five parts of wavelet power - the parts are representing the usual frequency
%          bands (delta, theta, beta, gamma).
%       2. Variance  vectors computed with cluster analysis.
%       3. Burst parameters
%       4. Discriminated unit ('vdisc').
%       5. A cell array containing the burst limits for every structural element (substate).
%       6. Cophenetic coefficient ('ccc' - see COPHENET for details).
%   All these outputs are stored in a mat file.
%
%   ICA_NEWAVE cuts waveletparts corresponding with burst localizations out. These 
%   burst-related waveletparts are linearly attached together ('linburstwave') and 
%   plotted as an image. The plots are saved in jpeg file format. ICA_NEWAVE calculates
%   average along time of the burst-related waveletparts and computes the angles of these
%   so called 'burstwave' vectors. These are saved in a text file.
%   Burst first spike location correlated values of the Hilbert transformated EEG (Hilbert
%   angles) are also saved in a text file.
%   'Burstwave' vectors are attached together and coefficient of variation (CV) along
%   'burst number' is calculated. These CV values are plotted against scale for every
%   substate determined by cluster analysis.
%
%   Output of EEG analysis part:
%       7. A cell array containing the CV (vs. scale) values for every substate ('scalecv').
%   This output is also stored in a mat file.
%
%   See also ICA, ICA_BETA, ICA_BETA2, ICA_BETA2B, ICA_WAVE, ICA_GUI2, ICA_GUI2B, 
%   WAVEPHASE_FOR_ICA, WAVELET, BURST_CLS, BURST_CLS_AUTO and CLUSTER.

% Input arguments check
error(nargchk(0,1,nargin));
if nargin == 0
    d = 20;
end

% Directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end
where = [DATAPATH,'DATA\analysenow3\'];    %Here are the data files
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'ICA\ica_newave\']);  %Here are the results
if ~isdir('Ward')   %create subdirectories
    mkdir Ward
end
if ~isdir('single')
    mkdir single
end
if ~isdir('average')
    mkdir average
end
if ~isdir('thres')
    mkdir thres
end
if ~isdir('wavevectors')
    mkdir wavevectors
end
if ~isdir('angles')
    mkdir angles
end
if ~isdir('wavelet_settings')
    mkdir wavelet_settings
end

% Text files
npt = input('Discard existing content or append data while writing ica_angles.txt? /discard:  ENTER, append: a/','s');
npt = input('Discard existing content or append data while writing ica_hilbert_angles.txt? /discard:  ENTER, append: a/','s');
if isempty(npt)
    cd angles
    fid_w = fopen('ica_angles.txt','w');
    fid_h = fopen('ica_hilbert_angles.txt','w');
    cd ..
elseif npt == 'a'
    cd angles
    fid_w = fopen('ica_angles.txt','a');
    fid_h = fopen('ica_hilbert_angles.txt','a');
    cd ..
else
    error('Unexpected answer for input.')
end

% Progress indicator
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
list = {};      %initialize list of analysed files
for o = 1:sf    %"CELL CYCLE"
    fname = files(o).name;
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
    seglen = 2000000;
    segno = floor(meret/seglen) + 1;
    for seg = 1:segno   %"REGISTRATION SEGMENT CYCLE"
        datinx1 = (seg - 1) * seglen + 1;      %first point of the interval
        datinx2 = min(seg*seglen,meret);       %last point of the interval
        b_imp(fname,where,data,datinx1,datinx2);
        
        islist = 1;     %initialize 'islist' for the decision of adding the cell name to the output list

% Loc. of first and last points of theta intervals
        if segno == 1
            fln = ['THETA_SEGMENTS_',filenam];
            ff = fullfile(DATAPATH,'Wavelet\segments',fln);
            load(ff)
        else
            fln = ['THETA_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'Wavelet\segments_long',fln);
            load(ff)
        end
        th_int = ThetaSegments;
        th_int1 = th_int(1,:) + 200;    %leaving 200 - 200 ms edges of each theta segment
        th_int2 = th_int(2,:) - 200;
        difth = th_int2 - th_int1;
        fd = find(difth<50000);         %leaving segments shorter than 5 sec.
        th_int1(fd) = [];
        th_int2(fd) = [];
        th_int = [th_int1; th_int2];
        th_int = th_int(:)';
        if isempty(th_int)      %FIRST CHECKPOINT: theta interval length criterium
            disp(['Problem with cell ' fn_noext ': no theta interval longer than 5 sec. was found.'])
            islist = 0;     %avoiding adding the cell name to the output cell list
            continue
        end
    
% Thresholding
        [T,H] = b_thres;
        kuszob = T;
        cd thres
        eval(['saveas(H,''THRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.fig'')']);
        cd ..
    
% Discrimination                
        b_disc(kuszob);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};
        kuszob = DISC{4};
        instfrek = DISC{5};
        isi = DISC{6};
    
% Selection of the corresponding unit segments
        [th_out_int,th_vd] = b_segmentselector(th_int,vdisc);
        if ~isempty(th_int) & isempty(th_out_int)   %SECOND CHECKPOINT: spike number criterium
            disp(['Problem with cell ' fn_noext ': no theta interval with more than 100 spikes was found.'])
            islist = 0;     %avoiding adding the cell name to the output cell list
            continue
        end
        pend = length(th_vd);
        vd = th_vd;
        out_int = th_out_int;
        st = 'theta';
    
% Wavelet transformation
        [wave,f,newstep] = b_waveletcall_for_applications;
        sw1 = size(wave,1);
        ScaleVector = f(1:sw1);
        Newstep = newstep;
        cd wavelet_settings
        save newstep Newstep
        save f ScaleVector
        cd ..
        
% Creating "wavelet vectors"
        fnd = find(f<1);
        pwind1 = fnd(1);
        fnd = find(f<3);
        pwind2 = fnd(1);
        fnd = find(f<6);
        pwind3 = fnd(1);
        fnd = find(f<20);
        pwind4 = fnd(1);
        fnd = find(f<50);
        pwind5 = fnd(1);
        
        lenw = size(wave,2);
        wavevec = zeros(4,lenw);                        %Wavevec contains the mean along scale of five parts of wavelet power
        wavevec(1,:) = mean(abs(wave(pwind5:pwind4-1,:)).^2);      %20 - 50 Hz
        wavevec(2,:) = mean(abs(wave(pwind4:pwind3-1,:)).^2);     %6 - 20 Hz
        wavevec(3,:) = mean(abs(wave(pwind3:pwind2-1,:)).^2);    %3 - 6 Hz
        wavevec(4,:) = mean(abs(wave(pwind2:pwind1-1,:)).^2);   %1 - 3 Hz
        
% Cluster analysis
        for p = 1:pend      %"THETA SEGMENT CYCLE"
            i_frst = out_int(2*p-1);    %create indeces
            i_scnd = out_int(2*p);
            ind1 = i_frst + (datinx1 - 1);
            ind2 = i_scnd + (datinx1 - 1);
            lgt = i_scnd - i_frst;
            preintstart = i_frst / newstep;
            preintend = i_scnd / newstep;
            intstart = floor(preintstart);
            intend = ceil(preintend);
            
            for lnk = 1:3       %"METHOD CYCLE"
                switch lnk
                case 1
                    method = 'ward';
                    cd Ward
                case 2
                    method = 'single';
                    cd single
                case 3
                    method = 'average';
                    cd average
                end
            
                [IntraBurstIvCv,ExtraBurstIvCv,InterBurstIvCv,FirstSpikeCv,AllFirstSpikeCv,...
                        Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,...
                        Burst,Ccc] = b_itclust(vd{p},d,method);     %cluster analysis and burst parameters
                
% Return plot & Vdisc
                ivs = diff(vd{p});
                livs = length(ivs);
                ReturnPlotXData = ivs(1:livs-1);
                ReturnPlotYData = ivs(2:livs);
                
                Vdisc = vd{p};
                
% Lomb periodogram
                fdi = find(diff(IntraBurstIvCv));    %find the different states
                fdi = fdi + 1;
                lenfdi = length(fdi);
%                 Pxx = cell(1,lenfdi);
%                 Pyy = cell(1,lenfdi);
%                 Z = cell(1,lenfdi);
%                 for c = 1:lenfdi
%                     dec = fdi(c);
%                     bas = [];
%                     for bno = 1:size(Burst{dec},2)
%                         b = Vdisc(Burst{dec}(1,bno):Burst{dec}(2,bno));
%                         bas = [bas b];
%                     end
%                     lenu = lgt;
%                     vdisc_for_lomb = Vdisc;
%                     [Pxx{c},Pyy{c},jmax,prob,Z{c},effm,period_handle] = b_lomb_period(vdisc_for_lomb,bas,lenu);
%                 end
                
% EEG analysis part
                ScaleCv = wave_segmentation(Burst,vd{p},wave,lenw,f,eeg,newstep,...
                    fdi,filenam,st,method,p,pend,ind1,ind2,fid_w,fid_h);
                            
% Saving
                str = ['THETA_ICA_UNIT_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' ReturnPlotXData ' 'ReturnPlotYData']);
                eval(['save ' str ' IntraBurstIvCv ' 'ExtraBurstIvCv ' 'InterBurstIvCv ' 'FirstSpikeCv ' 'AllFirstSpikeCv ' '-append']);
                eval(['save ' str ' Vdisc ' 'Burst ' 'Ccc ' '-append']);
                str = ['THETA_ICA_PARAM_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' Burstiness ' 'IntraBurstFrequency ' 'IntraBurstSpikeNumber ' 'BurstLength ' 'BurstFrequency']);
                str = ['THETA_ICA_EEG_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' ScaleCv']);
                str = ['THETA_ICA_LOMB_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
%                 eval(['save ' str ' Pxx' ' Pyy' ' Z']);
                cd ..
            end     %end of "method cycle"
                        
            cd wavevectors
            Wavevec  = wavevec(:,intstart:intend);
            str = ['WAVEVECTOR_',filenam,'_',num2str(ind1),'_',num2str(ind2)];
            eval(['save ' str ' Wavevec']);
            cd ..
        end     %end of "theta segment cycle"
    end     %end of "registration segment cycle"
    clear wave
    waitbar(o/sf)   %Progress indicator
    if islist    %generate list of analysed files
        list{end+1} = [filenam,'_',num2str(ind1),'_',num2str(ind2)];
    end
    
% Close    
    close all
end     %end of "cell cycle"
close(wb);   %Close progress indicator
fclose(fid_w);
fclose(fid_h);

% Save list of analysed files
list = list';
save list list
cd(mmm);

% -----------------------------------------------------------------------------------------------------------------
function ScaleCv = wave_segmentation(Burst,vdisc,wave,lenw,f,eeg,newstep,...
    fdi,filenam,st,method,p,pend,ind1,ind2,fid_w,fid_h)

% Hilbert transformation of the eeg
ahee = angle(hilbert(eeg));             
ahee_uw = unwrap(ahee);
ahee = ahee*(180/pi);

% Wave segmentation regarding the burst localizations
lenfdi = length(fdi);   %number of different states                 
vbfew = cell(1,lenfdi); %allocation
flag = cell(1,lenfdi);
burstwave = cell(1,lenfdi);
cv_bvec = cell(1,lenfdi);
lin_burstwave = cell(1,lenfdi);
linbw = cell(1,lenfdi);
linbvec = cell(1,lenfdi);                
bang = cell(1,lenfdi);
bang_uw = cell(1,lenfdi);
bang_bin = cell(1,lenfdi);
bang_count = cell(1,lenfdi);
bins = [-168:12:180];
clr = zeros(1,lenfdi);
legend_matrix = cell(1,lenfdi);
% cv_handle = figure;
for b = 1:lenfdi                        %first element of intraburstiv_cv is zeros
                                        %b indicates the different substates: "SUBSTATE CYCLE"
    vbfew{b} = vdisc(Burst{fdi(b)});    %vbfew contains the loc. of the first skipes of bursts in its first line
                                        %and loc. of the last spikes of bursts in its second line considering
                                        %only those steps where Burst matrix changes
    if size(vbfew{b}) == [1 2]
        vbfew{b} = vbfew{b}';
    end
    bang{b} = ahee(vbfew{b}(1,:));      %Hilbert angle at loc. of first spikes
    bang_uw{b} = ahee_uw(vbfew{b}(1,:));
    [bang_count{b} bang_bin{b}] = hist(bang{b},bins);
    burstlength = vbfew{b}(2,:) - vbfew{b}(1,:) + 1;    %prepare for plotting: 'burstlength', 'flag' and allocation
    flag{b} = min(burstlength,1000);
    lenvbfew = size(vbfew{b},2);
    bw = cell(1,lenvbfew);
    bvec = cell(1,lenvbfew);
    lin1 = cell(1,lenvbfew);
    lin2 = cell(1,lenvbfew);
    lin3 = [];
    lin4 = [];
    u141 = size(wave,1);
    for v = 1:lenvbfew      %v indicates the different bursts within the substate: "BURST CYCLE"
        lfsp = floor(vbfew{b}(1,v)/newstep);    %first x coordinate of burst (on wavelet)
        llsp = ceil(vbfew{b}(2,v)/newstep);     %last x coordinate of burst (on wavelet)
        ff = floor(flag{b}(v)/newstep);         %create left flag                   
        cf = ceil(flag{b}(v)/newstep);          %create right flag
        preindex1 = max(1,lfsp-ff);             %left-flagged burst, cut at first point of wavelet
        preindex2 = min(lenw,llsp+cf);          %right-flagged burst, cut at last point of wavelet
        index1 = preindex1;
        index2 = preindex2;
        bw{v} = wave(:,index1:index2);      %bw contains a three times (but max 2*1000) longer 
                                            %part of wave for each burst with the burst time
                                            %in the center
        bvec{v} = mean((abs(bw{v}).^2),2);  %bvec contains one vector for each burst with the
                                            %mean values of bw magnitude at each scale
        linbw{b} = [linbw{b} bw{v}];        %linbw contains the bw matrices linearly after 
                                            %each other at every b index
        linbvec{b} = [linbvec{b} bvec{v}];  %linbvec contains the bvec arrays linearly after 
                                            %each other in order to get a burst - scale matrix
                                            %at every b index
        slbw = size(linbw{b},2);
        pre1 = [slbw; slbw];
        pre2 = [1; u141];
        lin3 = [lin3 pre1];         %deliminator between bursts, y coordinates (on wavelet)
        lin4 = [lin4 pre2];         %deliminator between bursts, x coordinates (on wavelet)
        lbw = size(bw{v},2);
        lf1 = lfsp - preindex1;     %length of left flag (on wavelet)
        lf2 = preindex2 - llsp;     %length of right flag (on wavelet)
        lc = lbw - lf1 - lf2;       %length of burst (on wavelet)
        first = slbw - lf2 - lc;    %first x coord. of burst (on wavelet)
        last = slbw - lf2;          %last x coord. of burst (on wavelet)
        lin1{v} = [3 u141-1 u141-1 3 3];  %box for v. burst, y coordinates (on wavelet)
        lin2{v} = [first first last last first];    %box for v. burst, x coordinates (on wavelet)
    end     %end of "burst cycle"
    burstwave{b} = bw;  %burstwave contains a cell array at every b index with one matrix 
                        %per burst cut from wave
    lin_burstwave{b} = abs(linbw{b}).^2;    %lin_burstwave contains the power of linbw (linearly
                                            %attached bw matrices) for each substate
    cv_bvec{b} = std(linbvec{b},0,2) ./ (mean(linbvec{b},2)+eps);
                        %cv_bvec contains the coefficients of variation of burst - scale 
                        %matrices (linbvec) by burst for every substate
    tt1 = [filenam(1:3) filenam(5:6) ' ' st ' ' method ' subsegment: ' num2str(p) '/' num2str(pend)];
    
% Image
    levels = 2 .^ [-4:0.1:7];
    levels = [0:1:10];
    H = figure;
    imagesc(lin_burstwave{b});
    tt2 = [filenam(1:3) filenam(5:6) ' ' st ' ' method ' subsegment: ' num2str(p) '/' num2str(pend)...
            ' substate: ' num2str(b) '/' num2str(lenfdi)];
    title(tt2);
    b_rescaleaxis('Y',f)
    eval(['saveas(H,''THETA_ICA_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'_substate',num2str(b),'.jpg'')']);
                    
% Computing the mean angle
    lunder2 = lenvbfew * (lenvbfew - 1) / 2;
    angs = zeros(1,lunder2);
    next = 1;
    for out = 1:lenvbfew
        for in = out+1:lenvbfew
            angs(next) = b_ang(bvec{out},bvec{in});
            next = next + 1;
        end
    end
    if ~isempty(angs)
        meanrad = atan(sum(sin(angs))/sum(cos(angs)));
        meanangle = meanrad / pi * 180;
    else
        meanangle = 'X';
    end
    wr = [tt2 ' MEAN ANGLE: '];
    wr_h = [tt2 ' HILBERT ANGLE:'];
    fprintf(fid_w,'%s\n %E\n',wr,meanangle);
    fprintf(fid_h,'%s\n %E\n',wr_h,bang{b});
end     %end of "substate cycle"

% Output argument
ScaleCv = cv_bvec;