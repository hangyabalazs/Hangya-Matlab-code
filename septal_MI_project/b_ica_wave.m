function b_ica_wave(d)
%ICA_WAVE  Iterative Cluster Analysis and wavelet average vector calculation.
%   This function uses two directories - one for the data files and one for the results. 
%   You have to specify these directories in the program code.
%   ICA_WAVE(D) creates different number of clusters from 2 to D and computes intraburst 
%   interval, interburst interval, "inter-fist-spike interval", "all-inter-first-spike-
%   interval" (distances of first spikes including single spikes as well) and extraburst
%   interval length variance at all number of clusters.
%   ICA_WAVE, by itself, uses D = 20 default value.
%
%   ICA_WAVE, unlike ICA_BETA creates only matricies, not figures. You can display the figures
%   using ICA_GUI2B.
%
%   ICA_WAVE applies the following burst definition: intraburst intervals are the interspike
%   intervals containted by the interspike interval cluster with the shortest interval.
%
%   This function works on one predefinied interval (320000 - 920000) and loads previously
%   saved threshold files. Note that in case you modify the interval, you should repeat the
%   thresholding!
%
%   The program runs on theta and notheta segments in the predefinied interval. The selection
%   is based on WAVEPHASE. The segments have to be longer than 5 sec. and contain more than
%   100 spikes.
%
%   Wavelet spectrum is calculated using WAVELET function.
%
%   Outputs: 
%       1. Wavelet average magnitude vectors ('Wavevec') contain the mean along scale of 
%          five parts of wavelet power - the parts are representing the usual frequency
%          bands (delta, theta, beta, gamma).
%       2. Variance  vectors computed with cluster analysis.
%   All these outputs are mat files.
%
%   ICA_WAVE cuts waveletparts corresponding with burst localizations out and calculates
%   average along time of them. It computes the angles of these so called 'burstwave' 
%   vectors. These are saved in a text file.
%   Burst first spike location correlated values of the Hilbert transformated EEG (Hilbert
%   angles) are also saved in a text file.
%   'Burstwave' vectors are attached together and coefficient of variation (CV) along
%   'burst number' is calculated. These CV values are plotted against scale for every
%   substate determined by cluster analysis.
%
%   See also ICA, ICA_BETA, ICA_BETA2, ICA_BETA2B, ICA_GUI2, ICA_GUI2B, WAVEPHASE_FOR_ICA, 
%   WAVELET, BURST_CLS, BURST_CLS_AUTO and CLUSTER.

% Warning
disp(' ');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp('!!!           Note that in case you modify the interval,              !!!');
disp('!!!              you should repeat the thresholding!                  !!!');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp(' ');

% Input arguments check
error(nargchk(0,1,nargin));
switch nargin
case 0
    d = 21;
case 1
    d = d + 1;
end;

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'DATA\analysenow3\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'ICA\ica_wave_proba\']);  %Here are the results
npt = input('Discard existing content or append data while writing ica_angles.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid_w = fopen('ica_angles.txt','w');
    fid_h = fopen('ica_hilbert_angles','w');
elseif npt == 'a',
    fid_w = fopen('ica_angles.txt','a');
    fid_h = fopen('ica_hilbert_angles','a');
else
    error('Unexpected answer for input.')
end;

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    data = load(ffnm);
    meret=size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2)==1,
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
    lend = length(data);

% Downsample, new sampling freq. = 200 Hz
    datinx1 = 320000;
    datinx2 = 920000;
%     for segm = 1:2
%         switch segm
%         case 1
%             start = 0;
%             datinx1 = 1; %first point of the interval
%             datinx2 = 600000; %last point of the interval
%         case 2
%             start = 600000;
%             datinx1 = 600001; %first point of the interval
%             datinx2 = lend; %last point of the interval
%         end
        wholeeeg = data(datinx1:datinx2,1);
        lenwe = length(wholeeeg);
        newstep = 50;
        resamp = 10000/newstep;
        sst = wholeeeg(1:newstep:lenwe);
        
% Standardization
        variance = std(sst)^2;
        sst = (sst - mean(sst)) / sqrt(variance);
        
% Wavelet transformation
        dt = 1 / resamp;    %resample on 100 Hz
        time = [0:length(sst)-1] * dt + 0;
        n = length(sst);              %v
        pad = 1;      
        dj = 0.04;    %v
        s0 = 2 * dt;
        j1 = ((1 / dj) * log2(n/2)) * 2;
        j1 = ceil(j1);    %v
        j = (0:j1);               %v
        s = s0 .* 2 .^ (j * dj);                %v
        omega0 = 6;               %v
        c = 4 * pi / (omega0 + sqrt(2+omega0^2));               %v
        fperiod = c .* s;   %v
        f = 1 ./ fperiod;               %v
        lag1 = 0.72;  
        mother = 'Morlet';
        [wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
        
% Creating "wavelet vectors"
        interestwave = abs(wave(1:142,:));
        lenw = size(wave,2);
        wavevec = zeros(4,lenw);                        %Wavevec contains the mean along scale of five parts of wavelet power
        wavevec(1,:) = mean(interestwave(1:32,:));      %20 - 50 Hz
        wavevec(2,:) = mean(interestwave(33:76,:));     %6 - 20 Hz
        wavevec(3,:) = mean(interestwave(77:101,:));    %3 - 6 Hz
        wavevec(4,:) = mean(interestwave(102:142,:));   %1 - 3 Hz
        clear interestwave
        
% Computing the input variables
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
        IN{4} = where1;     %pathname
        IN{5} = datinx1;
        IN{6} = datinx2;
        IN{7} = time;
        IN{8} = unit;
        IN{9} = dt;
        IN{10} = meret;
        IN{11} = mintafr;
        IN{12} = xlimit;
        
% Loading the saved threshold value
        directory = [DATAPATH,'Data\megawavedisc3\'];
        str = [directory,'MEGAWAVE_',filenam,'.mat'];
        load(str)
        if isempty(kuszob)
            disp(['Problem with cell ',fname(1:6),': no kuszob.'])
            continue
        end
    
% Discrimination                
        b_disc(kuszob);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};
        kuszob = DISC{4};
        instfrek = DISC{5};
        isi = DISC{6};
        
% Loc. of first and last points of theta intervals and notheta intervals
        [th_int,no_int] = b_wavephase_for_ica(wave,f,newstep);
        if isempty(th_int)
            disp(['Problem with cell ' fn_noext ': no theta interval longer than 5 sec. was found.'])
        end
        
% Selection of the corresponding unit segments
        [th_out_int,th_vd] = b_segmentselector(th_int,vdisc,datinx1);
        [no_out_int,no_vd] = b_segmentselector(no_int,vdisc,datinx1);
        if ~isempty(th_int) & isempty(th_out_int)
            disp(['Problem with cell ' fn_noext ': no theta interval with more than 100 spikes was found.'])
        end
        
        for segment_type = 1:2
            switch segment_type
            case 1
                pend = length(th_vd);
                vd = th_vd;
                out_int = th_out_int;
                st = 'theta';
            case 2
                pend = length(no_vd);
                vd = no_vd;
                out_int = no_out_int;
                st = 'no theta';
            end
            for p = 1:pend
                
% Return plot
                ivs = diff(vd{p});
                livs = length(ivs);
                ReturnPlotXData = ivs(1:livs-1);
                ReturnPlotYData = ivs(2:livs);
                
% Ward's clustering
                dmtx = zeros(livs,livs);
                iivs = ivs';
                dmtx(1,:) = ivs;
                dmtx(:,1) = iivs;
                dist = pdist2(dmtx);
                links = linkage(dist,'ward');
                
                intraburstivvar_norm = zeros(1,d);  %allocating some matrices
                extraburstivvar_norm = zeros(1,d);
                interburstivvar_norm = zeros(1,d);
                firstspikevar_norm = zeros(1,d);
                allfirstspikevar_norm = zeros(1,d);
                
                Burst = cell(1,d);  %allocating Burst cell
                
                miniv = min(ivs);   %finding the cluster containing the smallest interval
                for dec = 3:d,
                    c = cluster(links,dec);
                    cell_clusters = cell(dec-1,1);
                    doub_length = zeros(dec-1,1);
                    for t = 1:dec-1,
                        cell_clusters{t} = find(c==t);
                        doub_length(t) = length(cell_clusters{t});
                        fnd = find(ivs(cell_clusters{t}) == miniv);
                        if isempty(fnd) == 0,
                            miniv_clus = t;
                        end;
                    end;
                    
                    ivss = ivs(cell_clusters{miniv_clus});  %creating extraburstiv
                    extraburstiv = ivs(find(ivs>max(ivss)));
                    
                    liv = find(ivs<=max(ivss)); %finding the bursts
                    liv1 = [-1 liv]; 
                    liv = [liv 0];
                    bliv = find(liv~=liv1+1);
                    blivv = bliv(2:end)-1; 
                    bliv(end) = [];
                    burst = [liv(bliv);liv(blivv)+1];     %1st and last spikes of bursts
                    fsp = vd{p}(burst(1,:));              %localisation of 1st spikes of bursts
                    lsp = vd{p}(burst(2,:));              %localisation of last spikes of bursts
                    diffburst = vd{p}(burst(2,:)) - vd{p}(burst(1,:));
                    mdb = max(diffburst);
                    
                    Burst{dec} = burst; %filling Burst cell
                    
                    ssi = vd{p};    %ssi is for "single spikes included": it will contain the first spikes of bursts
                                    %and the single spikes as well
                    intraburstiv = [];
                    for j = 1:size(burst,2),    %computing intraburstiv and allfirstspike
                        b = vd{p}(burst(1,j):burst(2,j));
                        intraburstiv = [intraburstiv diff(b)];                   
                        ssi(burst(1,j)+1:burst(2,j)) = 0;
                    end
                    
                    fsp2 = fsp(2:end);  %computing interburstiv
                    lsp2 = lsp(1:end-1);
                    interburstiv = fsp2 - lsp2;
                    
% Variance of the normalized intraburst interval length
                    if length(intraburstiv) ~= 0
                        intraburstivvar_norm(dec) = var(intraburstiv/mean(intraburstiv));  
                    else intraburstivvar_norm(dec) = NaN;
                    end
                    %intraburstivvar_norm2 = var(ivss/mean(ivss));
                    %this - much easier - way of computing intraburstivvar produced an unexpected error:
                    %the first intraburst interval was skipped
                    
% Variance of the normalized extraburst interval length
                    if length(extraburstiv) ~= 0,
                        extraburstivvar_norm(dec) = var(extraburstiv/mean(extraburstiv));
                    else extraburstivvar_norm(dec) = NaN;
                    end;
                    
% Variance of the normalized interburst interval length
                    if length(interburstiv) ~= 0,
                        interburstivvar_norm(dec) = var(interburstiv/mean(interburstiv));
                    else interburstivvar_norm(dec) = NaN;
                    end;
                    
% Variance of the normalized inter-1st-spike interval length
                    dfsp = diff(fsp);
                    if length(dfsp) ~= 0,
                        firstspikevar_norm(dec) = var(dfsp/mean(dfsp));
                    else firstspikevar_norm(dec) = NaN;
                    end;
                    
% Varinance of the normalized all-inter-1st-spike interval length
                    dfssi = diff(find(ssi));
                    allfirstspikevar_norm(dec) = var(dfssi/mean(dfssi));    %Note, that dfssi cannot be empty!
                end;                
                
% Hilbert transformation of the eeg
                ahee = angle(hilbert(eeg));             
                ahee_uw = unwrap(ahee);
                ahee = ahee*(180/pi);
                
% Wave segmentation regarding the burst localizations
                fdi = find(diff(intraburstivvar_norm)); %find the different states
                                                        %number of different states
                fdi = fdi + 1;
                lenfdi = length(fdi);                   
                vbfew = cell(1,lenfdi);
                flag = cell(1,lenfdi);
                burstwave = cell(1,lenfdi);
                cv_bvec = cell(1,lenfdi);
                linbw = cell(1,lenfdi);
                linbvec = cell(1,lenfdi);                
                bang = cell(1,lenfdi);
                bang_uw = cell(1,lenfdi);
                bang_bin = cell(1,lenfdi);
                bang_count = cell(1,lenfdi);
                bins = [-168:12:180];
                clr = zeros(1,lenfdi);
                legend_matrix = cell(1,lenfdi);
                cv_handle = figure;
                for b = 1:lenfdi                        %first two elements of intraburstivvar_norm are zeros
                                                        %b indicates the different substates
                    vbfew{b} = vd{p}(Burst{fdi(b)});    %vbfew contains the loc. of the first skipes of bursts in its first line
                                                        %and loc. of the last spikes of bursts in its second line considering
                                                        %only those steps where Burst matrix changes
                    if size(vbfew{b}) == [1 2]
                        vbfew{b} = vbfew{b}';
                    end
                    bang{b} = ahee(vbfew{b}(1,:));      %Hilbert angle at loc. of first spikes
                    bang_uw{b} = ahee_uw(vbfew{b}(1,:));
                    [bang_count{b} bang_bin{b}] = hist(bang{b},bins);
                    burstlength = vbfew{b}(2,:) - vbfew{b}(1,:) + 1;
                    flag{b} = min(burstlength,1000);
                    lenvbfew = size(vbfew{b},2);
                    bw = cell(1,lenvbfew);
                    bvec = cell(1,lenvbfew);
                    lin1 = cell(1,lenvbfew);
                    lin2 = cell(1,lenvbfew);
                    lin3 = [];
                    lin4 = [];
                    u141 = size(wave,1);
                    for v = 1:lenvbfew      %v indicates the different bursts within the substate
                        lfsp = floor(vbfew{b}(1,v)/newstep);
                        llsp = ceil(vbfew{b}(2,v)/newstep);
                        ff = floor(flag{b}(v)/newstep);                    
                        cf = ceil(flag{b}(v)/newstep);
                        preindex1 = max(1,lfsp-ff);
                        preindex2 = min(lenw,llsp+cf);
                        index1 = preindex1;
                        index2 = preindex2;
                        bw{v} = wave(:,index1:index2);      %bw contains a three times longer (but max 2*1000) 
                                                            %part of wave for each burst with the burst time
                                                            %in the center
                        bvec{v} = mean(abs(bw{v}),2);       %bvec contains one vector for each burst with the
                                                            %mean values of bw magnitude at each scale
                        linbw{b} = [linbw{b} bw{v}];        %linbw contains the bw matrices linearly after 
                                                            %each other at every b index
                        linbvec{b} = [linbvec{b} bvec{v}];  %linbvec contains the bvec arrays linearly after 
                                                            %each other in order to get a burst - scale matrix
                                                            %at every b index
                        slbw = size(linbw{1,b},2);
                        pre1 = [slbw; slbw];
                        pre2 = [1; u141];
                        lin3 = [lin3 pre1];
                        lin4 = [lin4 pre2];
                        lbw = size(bw{1,v},2);
                        lf1 = lfsp - preindex1;
                        lf2 = preindex2 - llsp;
                        lc = lbw - lf1 - lf2;
                        first = slbw - lf2 - lc;
                        last = slbw - lf2;
                        lin1{v} = [1.5 u141-0.5 u141-0.5 1.5 1.5];
                        lin2{v} = [first first last last first];
                    end
                    burstwave{b} = bw;  %burstwave contains a cell array at every b index with one matrix 
                                        %per burst cut from wave
                    cv_bvec{b} = std(linbvec{b},0,2) ./ (mean(linbvec{b},2)+eps);
                                        %cv_bvec contains the coefficients of variation of burst - scale 
                                        %matrices (linbvec) by burst for every substate
                    clr(b) = (b - 1) / (lenfdi - 1);
                    semilogx(f,cv_bvec{b},'Color',[0 1-clr(b) clr(b)]);
                    hold on
                    tt1 = [filenam(1:3) filenam(5:6) ' ' st ' subsegment: ' num2str(p) '/' num2str(pend)];
                    title(tt1);
                    x_lim = xlim;
                    y_lim = ylim;
                    axis([1 f(1) y_lim(1) y_lim(2)]);
                    legend_matrix{b} = ['substate: ' num2str(b) '/' num2str(lenfdi)];
                    
% Contour plot
                    levels=2.^[-4:0.1:7];
                    levels = [0:1:10];
%                     H = figure;
%                     C = contourc(abs(linbw{b}));
                    tt2 = [filenam(1:3) filenam(5:6) ' ' st ' subsegment: ' num2str(p) '/' num2str(pend) ' substate: '...
                            num2str(b) '/' num2str(lenfdi)];
%                     title(tt2);
%                     for vv = 1:lenvbfew
%                         line(lin2{vv},lin1{vv},'Color','k')
%                     end
%                     line(lin3,lin4,'Color','w')
%                     b_hidelines
%                     set(H,'KeyPressFcn','b_callhidelines')
                    
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
%                     disp(meanangle)
                    wr = [tt2 ' MEAN ANGLE: '];
                    wr_h = [tt2 'Hilbert_' num2str(b)];
                    fprintf(fid_w,'%s %E\n',wr,meanangle);
%                     plot(bang_uw{b}); hold on;
%                     figure; bar(bang_bin{b},bang_count{b});
%                     disp('Hilbert angles');
%                     disp(bang{b}');
                    fprintf(fid_h,'%s %E\n',wr_h,bang{b});
                end
                
                figure(cv_handle)
                legend(legend_matrix)     %creating the legend for the scale - CV plot
                
% Creating matrices to be saved
                ind1 = out_int(2*p-1);
                ind2 = out_int(2*p);
                lgt = ind2 - ind1;
                i_frst = (ind1 - 320000);
                i_scnd = (ind2 - 320000);
                preintstart = i_frst / newstep;
                preintend = i_scnd / newstep;
                intstart = floor(preintstart);
                intend = ceil(preintend);

                IntraBurstIvVar = intraburstivvar_norm; 
                ExtraBurstIvVar = extraburstivvar_norm;
                InterBurstIvVar = interburstivvar_norm;
                FirstSpikeVar = firstspikevar_norm;
                AllFirstSpikeVar = allfirstspikevar_norm;
                Vdisc = vd{p};
                Time = [1:lgt];
                Wavevec  = wavevec(:,intstart:intend);
                
% Saving
                if strcmp(st,'theta')
                    eval(['save(''THETA_ICA_',num2str(ind1),'_',num2str(ind2),'_',filenam,'.mat'',''ReturnPlotXData'',''ReturnPlotYData'')']);
                    eval(['save(''THETA_ICA_',num2str(ind1),'_',num2str(ind2),'_',filenam,...
                            '.mat'',''IntraBurstIvVar'',''ExtraBurstIvVar'',''InterBurstIvVar'',''FirstSpikeVar'',''AllFirstSpikeVar'',''Vdisc'',''Burst'',''Time'',''-append'')']);
                    str = ['WAVEVECTOR_',num2str(ind1),'_',num2str(ind2),'_',filenam];
                    eval(['save ' str ' ' 'Wavevec']);
                end
%                 close all
            end
            waitbar(o/sf)   %Progress indicator
        end    
%     end
end
close(wb);   %Close progress indicator
fclose(fid_w);
fclose(fid_h);
cd(mmm);