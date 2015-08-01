function b_ica_phase(d)
%ICA_NEWAVE2  Iterative Cluster Analysis and wavelet average vector calculation.
%   This function uses two directories - one for the data files and one for the results. 
%   You have to specify these directories in the program code.
%
%   ICA_NEWAVE2(D) creates different number of clusters from 2 to D and computes intraburst 
%   interval, interburst interval, "inter-fist-spike interval", "all-inter-first-spike-
%   interval" (distances of first spikes including single spikes as well) and extraburst
%   interval length CV (coefficient of variation: std(X) / mean(X)) at all number of clusters.
%   ICA_NEWAVE2, by itself, uses D = 20 default value.
%   Burst parameters are also computed: burstiness (number of intraburst spikes / number of
%   all spikes), average intraburst frequency (calculated for each burst and averaged),
%   average intraburst spike number, average burst length, frequency of burst first spikes
%   (called 'burst frequency').
%   Clustering methods 'Ward', 'single' and 'average' were used (see LINKAGE for details).
%
%   ICA_NEWAVE2, unlike ICA_BETA creates only matricies, not figures. You can display the figures
%   using ICA_GUI3.
%
%   ICA_NEWAVE2 applies the following burst definition: intraburst intervals are the interspike
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
%   ICA_NEWAVE2 cuts waveletparts corresponding with burst localizations out. These 
%   burst-related waveletparts are linearly attached together ('linburstwave') and 
%   plotted as an image. The plots are saved in jpeg file format. ICA_NEWAVE2 calculates
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
%   The difference between ICA_NEWAVE and ICA_NEWAVE2 is that ICA_NEWAVE2 is optimalized
%   for memory usage.
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
cd([DATAPATH,'ICA\ica_newave3\']);  %Here are the results
create_subdir;         %create subdirectories

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
    data = b_load_data(ffnm);
    
% Computing the input variables
    meret = size(data,1);
    seglen = 2000000;
    segno = floor(meret/seglen) + 1;
    for seg = 1:segno   %"REGISTRATION SEGMENT CYCLE"
        datinx1 = (seg - 1) * seglen + 1;      %first point of the interval
        datinx2 = min(seg*seglen,meret);       %last point of the interval
        b_imp(fname,where,data,datinx1,datinx2);
        clear time unit phase      %free memory
        
        islist = 1;     %initialize 'islist' for the decision of adding the cell name to the output list

% Loc. of first and last points of theta intervals
        if segno == 1
            fln = ['THETA_SEGMENTS_',filenam];
            ff = fullfile(DATAPATH,'Wavelet\theta_segments',fln);
            load(ff)
        else
            fln = ['THETA_SEGMENTS_LONG',num2str(seg),'_',filenam];
            ff = fullfile(DATAPATH,'Wavelet\theta_segments_long',fln);
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
    
% Thresholding & Discrimination
        [T,segmlen] = b_thres4;
        
        b_disc2(T,segmlen);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};
        
        cd thres
        str = ['THETA_ICA_THRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'];
        eval(['save ' str ' T segmlen vdisc']);
        cd ..
        
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
        sw2 = size(wave,2);
        
        clear DISC      %free memory
        clear global DISC
        clear global IN
        
%         save vartemp1 d where files sf mmm fid_w_ward fid_h_ward fid_w_average...
%             fid_h_average fid_w_single fid_h_single list islist o fname filenam...
%             fn_noext fln data meret seglen segno seg datinx1 datinx2 th_int...
%             out_int vdisc vd pend st sw1 sw2 f newstep DATAPATH
%         save vartemp2 wave
%         clear
%         load vartemp2
        
        phase = atan2(imag(wave), real(wave));
        
        clear wave
%         load vartemp1
        
%         pieceno = 10;
%         segm = fix(sw2/pieceno);
%         inx2 = 0;
%         next = 1;
%         while inx2 < sw2
%             inx1 = max(1,inx2+1);
%             inx2 = min(sw2,inx1+segm);
%             wavefrag = wave(:,inx1:inx2);
%             str = ['save temp' num2str(next) ' wavefrag'];
%             eval(str)
%             clear wavefrag
%             next = next + 1;
%         end
%         clear wave
%         power = [];
%         phase = [];
%         for wsn = 1:next-1
%             str = ['load temp' num2str(wsn)];
%             eval(str)
%             powerfrag = (abs(wavefrag)) .^2;
%             power = [power powerfrag];
%             clear powerfrag
%             phasefrag = angle(wavefrag);
%             phase = [phase phasefrag];
%             clear phasefrag
%             clear wavefrag
%         end
        
%         pieceno = 5;
%         segm = fix(sw2/pieceno);
%         power = [];
%         phase = [];
%         while ~isempty(wave)
%             index1 = 1;
%             index2 = min(segm,size(wave,2));
%             wavefrag = wave(:,index1:index2);
%             powerfrag = (abs(wavefrag)) .^ 2;
%             phasefrag = angle(wavefrag);
%             clear wavefrag
%             wave(:,index1:index2) = [];
%             power = [power powerfrag];
%             phase = [phase phasefrag];
%         end
        
        lenw = size(phase,2);
                
        ScaleVector = f(1:sw1);
        Newstep = newstep;
        cd wavelet_settings
        save newstep Newstep
        save f ScaleVector
        cd ..
        
% Saving whole wavelet
        H = figure;
        save_whole_wavelet(H,phase,filenam,datinx1,datinx2,f,newstep)
                
% ANALYSIS based on cluster analysis
        for p = 1:pend      %"THETA SEGMENT CYCLE"
            i_frst = out_int(2*p-1);    %create indices
            i_scnd = out_int(2*p);
            ind1 = i_frst + (datinx1 - 1);
            ind2 = i_scnd + (datinx1 - 1);
            lgt = i_scnd - i_frst;
            preintstart = i_frst / newstep;
            preintend = i_scnd / newstep;
            intstart = floor(preintstart);
            intend = ceil(preintend);
            
% Saving segment wavelet
            save_segmented_wavelet(H,phase,intstart,intend,filenam,ind1,ind2,f,newstep)
            
% Single point spike triggered wavelet
            spstw_phase = phase(:,round(vd{p}/newstep));
            
            edges = [-pi : pi/32 : pi];
            spstw_zsolt = [];
            for fr = 1:sw1
                hst = histc(spstw_phase(fr,:),edges);
                spstw_zsolt = [spstw_zsolt; hst(1:end-1)];
            end
            
            cd spstw
            imagesc(spstw_phase);     %phase
            tt2 = ['SPSTW ' filenam(1:3) ' ' filenam(5:6) ' ' num2str(ind1) ' ' num2str(ind2)];
            title(tt2);
            b_rescaleaxis('Y',f)
            cd phase
            eval(['saveas(H,''SPSTW_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.jpg'')']);
            
            imagesc(spstw_zsolt);     %Zsolt's method
            tt2 = ['SPSTW ZSOLT ' filenam(1:3) ' ' filenam(5:6) ' ' num2str(ind1) ' ' num2str(ind2)];
            title(tt2);
            x_lim = get(gca,'XLim');
            y_lim = get(gca,'YLim');
            xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
            ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
            text(xcoord,ycoord,num2str(length(vd{p})),'Color','black','FontWeight','bold','FontSize',18);
            b_rescaleaxis('Y',f)
            xx = linspace(-pi,pi,length(edges)-1);
            b_rescaleaxis('X',xx)
            eval(['saveas(H,''SPSTW_ZSOLT_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.jpg'')']);
            cd ..
            cd ..
            
            clear spstw_phase aspstw_phase sspstw_phase   % free memory
                    
            % Clustering
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
            
                [IntraBurstIvCv,ExtraBurstIvCv,InterBurstIvCv,FirstSpikeCv,AllFirstSpikeCv,BurstLengthCv,...
                        Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,...
                        Burst,Ccc] = b_itclust2(vd{p},d,method);     %cluster analysis and burst parameters
                
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
                [ScaleCvPhase StwPhaseMean StwPhaseStd] =...
                    wave_segmentation(Burst,Vdisc,phase,lenw,f,newstep,fdi,filenam,st,method,p,pend,...
                    ind1,ind2,H,seg);
                            
% Saving
                str = ['THETA_ICA_UNIT_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' ReturnPlotXData ReturnPlotYData']);
                eval(['save ' str ' IntraBurstIvCv ExtraBurstIvCv InterBurstIvCv FirstSpikeCv AllFirstSpikeCv BurstLengthCv -append']);
                eval(['save ' str ' Vdisc Burst Ccc -append']);
                str = ['THETA_ICA_PARAM_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' Burstiness IntraBurstFrequency IntraBurstSpikeNumber BurstLength BurstFrequency']);
                str = ['THETA_ICA_EEG_PHASE_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' ScaleCvPhase']);
                str = ['THETA_ICA_STW_PHASE_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
                eval(['save ' str ' StwPhaseMean StwPhaseStd']);
%                 str = ['THETA_ICA_LOMB_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.mat'];
%                 eval(['save ' str ' Pxx' ' Pyy' ' Z']);
                cd .. ;
            end     %end of "method cycle"
                        
            if islist    %generate list of analysed files
                list{end+1} = [filenam,'_',num2str(ind1),'_',num2str(ind2)];
            end
        end     %end of "theta segment cycle"
    end     %end of "registration segment cycle"
    waitbar(o/sf)   %Progress indicator
    
% Close
    close all
end     %end of "cell cycle"
close(wb);   %Close progress indicator

% Save list of analysed files
list = list';
save phase_list list
cd(mmm);

% -----------------------------------------------------------------------------------------------------------------
function [ScaleCvPhase, StwPhaseMean, StwPhaseStd] =...
    wave_segmentation(Burst,Vdisc,phase,lenw,f,newstep,fdi,filenam,st,method,p,pend,...
    ind1,ind2,H,seg)

% Wave segmentation regarding the burst localizations
lenfdi = length(fdi);   %number of different states
cv_bvec_phase = cell(1,lenfdi);
bins = [-168:12:180];
clr = zeros(1,lenfdi);
legend_matrix = cell(1,lenfdi);
for b = 1:lenfdi                        %first element of intraburstiv_cv is zeros
                                        %b indicates the different substates: "SUBSTATE CYCLE"
    vbfew = Vdisc(Burst{fdi(b)});    %vbfew contains the loc. of the first spikes of bursts in its first line
                                     %and loc. of the last spikes of bursts in its second line considering
                                     %only those steps where Burst matrix changes
    if size(vbfew) == [1 2]
        vbfew = vbfew';
    end
    burstlength = vbfew(2,:) - vbfew(1,:) + 1;    %prepare for plotting: 'burstlength', 'flag' and allocation
    flag = min(burstlength,1000);
    lenvbfew = size(vbfew,2);
    bvec_phase = cell(1,lenvbfew);
    lin3 = [];
    lin4 = [];
    linbw_phase = [];
    linbvec_phase = [];
    u141 = size(phase,1);
    for v = 1:lenvbfew      %v indicates the different bursts within the substate: "BURST CYCLE"
        lfsp = floor(vbfew(1,v)/newstep);    %first x coordinate of burst (on wavelet)
        llsp = ceil(vbfew(2,v)/newstep);     %last x coordinate of burst (on wavelet)
        ff = floor(flag(v)/newstep);         %create left flag
        cf = ceil(flag(v)/newstep);          %create right flag
        preindex1 = max(1,lfsp-ff);             %left-flagged burst, cut at first point of wavelet
        preindex2 = min(lenw,llsp+cf);          %right-flagged burst, cut at last point of wavelet
        index1 = preindex1;
        index2 = preindex2;
        bw_phase = phase(:,index1:index2);   %bw_phase contains a three times (but max 2*1000) longer
                                             %part of wave phase for each burst with the burst time
                                             %in the center
        bvec_phase{v} = b_circular_mean(bw_phase,2);    %bvec_phase contains one vector for each burst with the
                                                        %mean values of bw_phase at each scale
        linbw_phase = [linbw_phase bw_phase];        %linbw_phase contains the bw_phase matrices linearly
                                                     %after each other
        linbvec_phase = [linbvec_phase bvec_phase{v}];  %linbvec_phase contains the bvec_phase arrays
                                                        %linearly after each other in order to get a
                                                        %burst - scale phase matrix
        slbw = size(linbw_phase,2);
        pre1 = [slbw; slbw];
        pre2 = [1; u141];
        lin3 = [lin3 pre1];         %deliminator between bursts, y coordinates (on wavelet)
        lin4 = [lin4 pre2];         %deliminator between bursts, x coordinates (on wavelet)
        lbw = size(bw_phase,2);
        lf1 = lfsp - preindex1;     %length of left flag (on wavelet)
        lf2 = preindex2 - llsp;     %length of right flag (on wavelet)
        lc = lbw - lf1 - lf2;       %length of burst (on wavelet)
        first = slbw - lf2 - lc;    %first x coord. of burst (on wavelet)
        last = slbw - lf2;          %last x coord. of burst (on wavelet)
        lin1 = [3 u141-1 u141-1 3 3];  %box for v. burst, y coordinates (on wavelet)
        lin2 = [first first last last first];    %box for v. burst, x coordinates (on wavelet)
    end     %end of "burst cycle"
    cv_bvec_phase{b} = b_circular_std(linbvec_phase,2) ./ (b_circular_mean(linbvec_phase,2)+eps);
                                    %cv_bvec_phase contains the coefficients of variation of burst - scale
                                    %phase matrices (linbvec_phase) by burst for every substate
    tt1 = [filenam(1:3) filenam(5:6) ' ' st ' ' method ' subsegment: ' num2str(p) '/' num2str(pend)];

% Image
    levels = 2 .^ [-4:0.1:7];
    levels = [0:1:10];

    imagesc(linbw_phase);     %phase
    tt2 = [filenam(1:3) filenam(5:6) ' ' st ' ' method ' subsegment: ' num2str(p) '/' num2str(pend)...
        ' substate: ' num2str(b) '/' num2str(lenfdi)];
    title(tt2);
    b_rescaleaxis('Y',f)
    cd phase
    eval(['saveas(H,''BURSTWAVE_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'_substate',num2str(b),'.jpg'')']);
    cd ..

% Spike trigerred wavelet
    dec = fdi(b);
    bas = [];
    for bno = 1:size(Burst{dec},2)
        bs = Vdisc(Burst{dec}(1,bno):Burst{dec}(2,bno));
        bas = [bas bs];
    end

    stw_phase = phase(:,round(bas/newstep));
    astw_phase{b} = b_circular_mean(stw_phase,2);
    sstw_phase{b} = b_circular_std(stw_phase,2);

    cd stw
    imagesc(stw_phase);     %phase
    tt2 = [filenam(1:3) filenam(5:6) ' ' st ' ' method ' subsegment: ' num2str(p) '/' num2str(pend)...
        ' substate: ' num2str(b) '/' num2str(lenfdi)];
    title(tt2);
    b_rescaleaxis('Y',f)
    cd phase
    eval(['saveas(H,''STW_',method,'_',filenam,'_',num2str(ind1),'_',num2str(ind2),'_substate',num2str(b),'.jpg'')']);
    cd ..
    cd ..
end

% Output
ScaleCvPhase = cv_bvec_phase;
StwPhaseMean = astw_phase;
StwPhaseStd = sstw_phase;

% ----------------------------------------------------------------------------------
function save_whole_wavelet(H,phase,filenam,datinx1,datinx2,f,newstep)

% Plot phase
cd wavelet
imagesc(phase);
stt = 'WAVELET WHOLE PHASE';
tt2 = [filenam(1:3) filenam(5:6) ' ' num2str(datinx1) ' ' num2str(datinx2) ' ' stt];
title(tt2);
lenp = datinx2 - datinx1 + 1;
sw2 = size(phase,2);
wavetime = (([1:sw2] - 1) * newstep + 1) / 10000;
ff = round(f*100) / 100;
time = round(wavetime*100) / 100;
b_rescaleaxis('Y',ff)
b_rescaleaxis('X',time)
cd phase
eval(['saveas(H,''WAVELET_PHASE_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.jpg'')']);
cd ..
cd ..

% ----------------------------------------------------------------------------------
function save_segmented_wavelet(H,phase,intstart,intend,filenam,ind1,ind2,f,newstep)

% Plot phase
cd wavelet
imagesc(phase(:,intstart:intend));
stt = 'WAVELET PHASE';
lenp = intend - intstart;
wavetime = [0:lenp-1] * newstep / 10000;
ff = round(f*100) / 100;
time = round(wavetime*100) / 100;
b_rescaleaxis('Y',ff)
b_rescaleaxis('X',time)
cd phase
eval(['saveas(H,''WAVELET_PHASE_',filenam,'_',num2str(ind1),'_',num2str(ind2),'.jpg'')']);
cd ..
cd ..

% ----------------------------------------------------------------------------------
function create_subdir

% Create subdirectories
if ~b_isdir2('Ward')
    mkdir Ward
end
cd Ward
if ~b_isdir2('phase')
    mkdir phase
end
if ~b_isdir2('stw')
    mkdir stw
end
cd stw
if ~b_isdir2('phase')
    mkdir phase
end
cd ..
cd ..

if ~b_isdir2('single')
    mkdir single
end
cd single
if ~b_isdir2('phase')
    mkdir phase
end
if ~b_isdir2('stw')
    mkdir stw
end
cd stw
if ~b_isdir2('phase')
    mkdir phase
end
cd ..
cd ..

if ~b_isdir2('average')
    mkdir average
end
cd average
if ~b_isdir2('phase')
    mkdir phase
end
if ~b_isdir2('stw')
    mkdir stw
end
cd stw
if ~b_isdir2('phase')
    mkdir phase
end
cd ..
cd ..

if ~b_isdir2('wavelet')
    mkdir wavelet
end
cd wavelet
if ~b_isdir2('phase')
    mkdir phase
end
cd ..

if ~b_isdir2('spstw')
    mkdir spstw
end
cd spstw
if ~b_isdir2('phase')
    mkdir phase
end
cd ..

if ~b_isdir2('thres')
    mkdir thres
end
if ~b_isdir2('wavelet_settings')
    mkdir wavelet_settings
end