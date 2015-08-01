function b_ica_beta2(d)
%ICA_BETA2  Beta2 version of Iterative Cluster Analysis.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code.
%   ICA_BETA2(D) creates different number of clusters from 2 to D and computes intraburst 
%   interval, interburst interval, "inter-fist-spike interval", "all-inter-first-spike-
%   interval" (distances of first spikes including single spikes as well) and extraburst
%   interval length variance at all number of clusters.
%   ICA_BETA2, by itself, uses D = 20 default value.
%
%   ICA_BETA2, unlike ICA_BETA creates only matricies, not figures. You can display the figures
%   using ICA_GUI2.
%
%   ICA_BETA2 applies the following burst definition: intraburst intervals are the interspike
%   intervals containted by the interspike interval cluster with the shortest interval.
%
%   This function works on the predefnied evoked theta, spontanous theta and no theta intervals.
%   See ICA_BETA2B or ICA_WAVE for other solutions.
%
%   See also ICA, ICA_BETA, ICA_BETA2B, ICA_WAVE, ICA_GUI2, ICA_GUI2B, BURST_CLS, BURST_CLS_AUTO
%   and CLUSTER.

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
where1 = [DATAPATH,'DATA\analysenow2\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'DATA\analysenow1\'];    %Here are the datinx 
mmm = pwd;
cd([DATAPATH,'ICA\ica_beta_2004\']);  %Here are the results
% npt = input('Discard existing content or append data while writing ica_beta.txt? /discard:  ENTER, append: a/','s');
% if isempty(npt),
%     fid = fopen('ica_beta.txt','w');
% elseif npt == 'a',
%     fid = fopen('ica_beta.txt','a');
% else
%     error('Unexpected answer for input.')
% end;

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    meret=size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2)==1,
        data=[data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
    
% Searching for the matching datinx file
    filenam = fname(1:6);
    filename = [filenam,'.txt'];
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
    if isempty(v)
        disp(['Problem with cell ',fname(1:6),': no datinx.'])
        continue
    end
    ww = 0;
    
% Computing the input variables
    for p = 1:3,
        datinx1 = v(2*p-1); %first point of the interval
        datinx2 = v(2*p); %last point of the interval
        
        dt = 0.0001;
        mintafr = 1 / dt;
        if datinx1(1) ~= 0 & datinx2(1) ~= 0,
            ww = ww + 1;
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
            IN{4} = where1; %pathname
            IN{5} = datinx1;
            IN{6} = datinx2;
            IN{7} = time;
            IN{8} = unit;
            IN{9} = dt;
            IN{10} = meret;
            IN{11} = mintafr;
            IN{12} = xlimit;
            
% Discrimination                
            kuszob = v(6+ww);
            try
                b_disc(kuszob);
            catch
                str = ['Error during discrimination of cell ',fn_noext,'.'];
                disp(str);
                break
            end
            global DISC
            id = DISC{1};
            output = DISC{2};
            vdisc = DISC{3};
            kuszob = DISC{4};
            instfrek = DISC{5};
            isi = DISC{6};
            
% Spike number criterium            
            if length(vdisc) > 10 & length(vdisc) > d,
                
% Return plot
                ivs = diff(vdisc);
                livs = length(ivs);
                ReturnPlotXData = ivs(1:livs-1);
                ReturnPlotYData = ivs(2:livs);
                switch p
                case 1
                    eval(['save(''NO_THETA_ICA_',filenam,'.mat'',''ReturnPlotXData'',''ReturnPlotYData'')']);
                case 2
                    eval(['save(''SP_THETA_ICA_',filenam,'.mat'',''ReturnPlotXData'',''ReturnPlotYData'')']);
                case 3
                    eval(['save(''EV_THETA_ICA_',filenam,'.mat'',''ReturnPlotXData'',''ReturnPlotYData'')']);
                end;
                
% Ward's clustering
                dmtx = zeros(livs,livs);
                iivs = ivs';
                dmtx(1,:) = ivs;
                dmtx(:,1) = iivs;
                dist = pdist2(dmtx);
                links = linkage(dist,'ward');
%                 [q w] = dendrogram(links);
                
                intraburstivvar_norm = zeros(1,d);  %allocating some matrices
                extraburstivvar_norm = zeros(1,d);
                interburstivvar_norm = zeros(1,d);
                firstspikevar_norm = zeros(1,d);
                allfirstspikevar_norm = zeros(1,d);
%                 intraburstivvar = zeros(1,d);
%                 extraburstivvar = zeros(1,d);
%                 interburstivvar = zeros(1,d);
%                 firstspikevar = zeros(1,d);

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
                    fsp = vdisc(burst(1,:));              %localisation of 1st spikes of bursts
                    lsp = vdisc(burst(2,:));              %localisation of last spikes of bursts
                    diffburst = vdisc(burst(2,:)) - vdisc(burst(1,:));
                    mdb = max(diffburst);
                    
                    Burst{dec} = burst; %filling Burst cell
                    
                    LE = zeros(1,length(burst));    %LE is for intraburstiv computation
                    ssi = vdisc;    %ssi is for "single spikes included": it will contain the loc. of the first 
                                    %spikes of bursts and the single spikes as well
                    for j = 1:size(burst,2),    %computing intraburstiv and allfirstspike
                        b = vdisc(burst(1,j):burst(2,j));
                        eval(['B',num2str(j),'=diff(b);']);
                        eval(['LE(j)=length(B',num2str(j),');']);
                        
                        ssi(burst(1,j)+1:burst(2,j)) = 0;
                    end;
                    intraburstiv = zeros(1,sum(LE));
                    next = 1;
                    for j = 1:length(burst),
                        for k = 1:LE(j),
                            eval(['intraburstiv(next) = B',num2str(j),'(k);']);
                            next = next + 1;
                        end;
                    end;
                    
                    fsp2 = fsp(2:end);  %computing interburstiv
                    lsp2 = lsp(1:end-1);
                    interburstiv = fsp2 - lsp2;
                    
% Variance of the normalized intraburst interval length
                    if length(intraburstiv) ~= 0,
                        intraburstivvar_norm(dec) = var(intraburstiv/mean(intraburstiv));  
                    else intraburstivvar_norm(dec) = NaN;
                    end;
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
                    
% Variance of the normalized all-inter-1st-spike interval length
                    dfssi = diff(find(ssi));
                    allfirstspikevar_norm(dec) = var(dfssi/mean(dfssi));    %Note, that dfssi cannot be empty
                    
%                 fssi=ssi(find(ssi));
%                 tfssi=zeros(1,length(unit));
%                 tfssi(fssi)=1;
%                 figure
%                 subplot(2,1,1);
%                 plot(time,unit)
%                 subplot(2,1,2);
%                 plot(time,tfssi)
%                 str = [dec];
%                 title(str)
                    
                end;
                
% Creating matrices to be saved                                
                IntraBurstIvVar = intraburstivvar_norm; 
                ExtraBurstIvVar = extraburstivvar_norm;
                InterBurstIvVar = interburstivvar_norm;
                FirstSpikeVar = firstspikevar_norm;
                AllFirstSpikeVar = allfirstspikevar_norm;
                Vdisc = vdisc;
                Time = time;
                
% Saving
                switch p    
                case 1
                    eval(['save(''NO_THETA_ICA_',filenam,...
                            '.mat'',''IntraBurstIvVar'',''ExtraBurstIvVar'',''InterBurstIvVar'',''FirstSpikeVar'',''AllFirstSpikeVar'',''Vdisc'',''Burst'',''Time'',''-append'')']);
                case 2
                    eval(['save(''SP_THETA_ICA_',filenam,...
                            '.mat'',''IntraBurstIvVar'',''ExtraBurstIvVar'',''InterBurstIvVar'',''FirstSpikeVar'',''AllFirstSpikeVar'',''Vdisc'',''Burst'',''Time'',''-append'')']);
                case 3
                    eval(['save(''EV_THETA_ICA_',filenam,...
                            '.mat'',''IntraBurstIvVar'',''ExtraBurstIvVar'',''InterBurstIvVar'',''FirstSpikeVar'',''AllFirstSpikeVar'',''Vdisc'',''Burst'',''Time'',''-append'')']);
                end;
            end;
        end;
    end;
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
% fclose(fid);
cd(mmm);