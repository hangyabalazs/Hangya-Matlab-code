function b_ica_beta(d)
%ICA_BETA  Beta version of iterative cluster analysis.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code.
%   ICA_BETA(D) creates different number of clusters from 2 to D and plots intraburst interval,
%   interburst interval, "inter-fist-spike interval" and extraburst interval length variance 
%   versus number of clusters. It also saves the variance values in a text file.
%   ICA_BETA, by itself, uses D = 20 default value.
%
%   ICA_BETA applies the following burst definition: intraburst intervals are the interspike
%   intervals containted by the interspike interval cluster with the shortest interval.
%
%   See also ICA, BURST_CLS, BURST_CLS_AUTO and CLUSTER.

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
where1 = [DATAPATH,'Data\analysenow3\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'Data\analysenow1\'];    %Here are the datinx 
mmm = pwd;
cd([DATAPATH,'ICA\ica_beta_figs\']);  %Here are the results
npt = input('Discard existing content or append data while writing ica_beta.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('ica_beta.txt','w');
elseif npt == 'a',
    fid = fopen('ica_beta.txt','a');
else
    error('Unexpected answer for input.')
end;

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
    cl = fname(1:6);
    filename = [cl,'.txt'];
    pont = findstr(filename,'.');
    filenam = filename(1:pont(1)-1);
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
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
            IN{4} = pathname;
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
            b_disc(kuszob);
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
                d1 = ivs(1:livs-1);
                d2 = ivs(2:livs);
                h1 = figure;
                plot(d1,d2,'.')
                title('return plot of ISIs');
                xlabel('nth ISI');
                ylabel('(n+1)th ISI');
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' RETURN PLOT']);
                    eval(['saveas(h1,''NO_THETA_RETURN_PLOT_',filenam,'.fig'')']);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' RETURN PLOT']);
                    eval(['saveas(h1,''SP_THETA_RETURN_PLOT_',filenam,'.fig'')']);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' RETURN PLOT']);
                    eval(['saveas(h1,''EV_THETA_RETURN_PLOT_',filenam,'.fig'')']);
                end;
                % k=waitforbuttonpress;
                % if k==0 close; end;
                
% Ward's clustering
                dmtx = zeros(livs,livs);
                iivs = ivs';
                dmtx(1,:) = ivs;
                dmtx(:,1) = iivs;
                dist = pdist2(dmtx);
                links = linkage(dist,'ward');
                [q w] = dendrogram(links);
                
                intraburstivvar_norm = zeros(1,d);
                extraburstivvar_norm = zeros(1,d);
                interburstivvar_norm = zeros(1,d);
                firstspikevar_norm = zeros(1,d);
                intraburstivvar = zeros(1,d);
                extraburstivvar = zeros(1,d);
                interburstivvar = zeros(1,d);
                firstspikevar = zeros(1,d);
                miniv = min(ivs);
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
%                     max_length = find(doub_length==max(doub_length));
%                                         
%                     
%                     binum=fix(0.606+0.4*(livs-1));
%                     
%                     if length(max_length) > 1,
%                         mm = zeros(1,length(max_length));
%                         for z = 1:length(max_length),
%                             mm(z) = min(ivs(cell_clusters{max_length(z)}));
%                         end;
%                         fmm = find(mm==min(mm));
%                         max_length = max_length(fmm);
%                     end;
%                     
%                     ivss = ivs(cell_clusters{max_length});

                    ivss = ivs(cell_clusters{miniv_clus});
                    extraburstiv = ivs(find(ivs>max(ivss)));
                    
                    
                    liv = find(ivs<=max(ivss));
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
                    
                    dt = 0.0001;
                    time = [0:length(unit)-1]*dt;
                    close all
                    h = figure;
                    plot(time,unit,'m');
                    hold on;
                    mi = kuszob; 
                    ma = kuszob + 0.5;
                    LE = zeros(1,length(burst));
                    for j = 1:size(burst,2),
                        rajz = [time(vdisc(burst(1,j))) time(vdisc(burst(2,j)));mi ma];
                        line(rajz(1,:),rajz(2,:),'color','b');
                        b = vdisc(burst(1,j):burst(2,j));
                        eval(['B',num2str(j),'=diff(b);']);
                        eval(['LE(j)=length(B',num2str(j),');']);
                    end;
                    intraburstiv = zeros(1,sum(LE));
                    next = 1;
                    for j = 1:length(burst),
                        for k = 1:LE(j),
                            eval(['intraburstiv(next) = B',num2str(j),'(k);']);
                            next = next + 1;
                        end;
                    end;
                    
                    fsp2 = fsp(2:end);
                    lsp2 = lsp(1:end-1);
                    interburstiv = fsp2 - lsp2;
                    
% Variance of the normalized intraburst interval length
                    if length(intraburstiv) ~= 0,
                        intraburstivvar_norm(dec) = var(intraburstiv/mean(intraburstiv));  
                    else intraburstivvar_norm(dec) = NaN;
                    end;
                    %intraburstivvar_norm2 = var(ivss/mean(ivss));
                    
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
                    
% Variance of the intraburst interval length
                    intraburstivvar(dec) = var(intraburstiv);    
                    %intraburstivvar2 = var(ivss);
                    
% Variance of the extraburst interval length
                    extraburstivvar(dec) = var(extraburstiv);
                    
% Variance of the interburst interval length
                    interburstivvar(dec) = var(interburstiv);
                    
% Variance of the inter-1st-spike interval length
                    dfsp = diff(fsp);
                    firstspikevar(dec) = var(dfsp);
                     
% Variance of the normalized inter-1st-spike interval length
                    if length(dfsp) ~= 0,
                        firstspikevar_norm(dec) = var(dfsp/mean(dfsp));
                    else firstspikevar_norm(dec) = NaN;
                    end;

% Plotting and saving the variance values in a text file
                    hold off;
                    title('Bursts');
                    ymin = min(unit) - 0.5;
                    ymax = max(unit) + 2;
                    x_lim = xlim;
                    axis([x_lim(1) x_lim(2) ymin ymax]);
                    y_lim = ylim;
                    text(x_lim(1)+1,y_lim(2)-0.1,['{\itVariance of the intraburst interval length }',num2str(intraburstivvar(dec))]);
                    %                     text(x_lim(1)+1,y_lim(2)-0.3,['{\itVariance of the intraburst interval length2 }',num2str(intraburstivvar2)]);
                    text(x_lim(1)+1,y_lim(2)-0.3,['{\itVariance of the extraburst interval length }',num2str(extraburstivvar(dec))]);
                    text(x_lim(1)+1,y_lim(2)-0.5,['{\itVariance of the interburst interval length }',num2str(interburstivvar(dec))]);
                    text(x_lim(1)+1,y_lim(2)-0.7,['{\itVariance of the inter-1st-spike interval length }',num2str(firstspikevar(dec))]);
                    text(x_lim(1)+1,y_lim(2)-0.9,['{\itVariance of the normalized intraburst interval length }',num2str(intraburstivvar_norm(dec))]);
                    %                     text(x_lim(1)+1,y_lim(2)-0.9,['{\itVariance of the normalized intraburst interval length2 }',num2str(intraburstivvar_norm2)]);
                    text(x_lim(1)+1,y_lim(2)-1.1,['{\itVariance of the normalized extraburst interval length }',num2str(extraburstivvar_norm(dec))]);
                    text(x_lim(1)+1,y_lim(2)-1.3,['{\itVariance of the normalized interburst interval length }',num2str(interburstivvar_norm(dec))]);
                    text(x_lim(1)+1,y_lim(2)-1.5,['{\itVariance of the normalized inter-1st-spike interval length }',num2str(firstspikevar_norm(dec))]);
                    if length(vdisc) < 100,
                        text(x_lim(1)+1,y_lim(2)-1.8,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','r');
                    else
                        text(x_lim(1)+1,y_lim(2)-1.8,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','k');
                    end;
                    
                    pont = findstr(filename,'.');
                    filenam = filename(1:pont(1)-1);
                    switch p
                    case 1
                        title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' ICA {\itdec} = ',int2str(dec)]);
                        eval(['saveas(h,''NO_THETA_ICA_',filenam,'_dec',int2str(dec),'.fig'')']);
                        ffilenam = [filenam,' NO THETA'];
                        fprf = [ffilenam,' ', num2str(dec)];
                        fprintf(fid,'%s %s\n',fprf);
                        fprintf(fid,'\n',fprf);
                        fprf = [num2str(intraburstivvar(dec)),' ',num2str(extraburstivvar(dec)),' ',num2str(interburstivvar(dec)),' ',...
                                num2str(firstspikevar(dec)),' ',num2str(intraburstivvar_norm(dec)),' ',num2str(extraburstivvar_norm(dec)),...
                                ' ',num2str(interburstivvar_norm(dec)),' ',num2str(firstspikevar_norm(dec))];
                        fprintf(fid,'%s %s %s %s %s %s %s %s\n',fprf);
                        fprintf(fid,'\n',fprf);
                    case 2
                        title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' ICA {\itdec} = ',int2str(dec)]);
                        eval(['saveas(h,''SP_THETA_ICA_',filenam,'_dec',int2str(dec),'.fig'')']);
                        ffilenam = [filenam,' SP THETA'];
                        fprf = [ffilenam,' ', num2str(dec)];
                        fprintf(fid,'%s %s\n',fprf);
                        fprintf(fid,'\n',fprf);
                        fprf = [num2str(intraburstivvar(dec)),' ',num2str(extraburstivvar(dec)),' ',num2str(interburstivvar(dec)),' ',...
                                num2str(firstspikevar(dec)),' ',num2str(intraburstivvar_norm(dec)),' ',num2str(extraburstivvar_norm(dec)),...
                                ' ',num2str(interburstivvar_norm(dec)),' ',num2str(firstspikevar_norm(dec))];
                        fprintf(fid,'%s %s %s %s %s %s %s %s\n',fprf);
                        fprintf(fid,'\n',fprf);
                    case 3
                        title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' ICA {\itdec} = ',int2str(dec)]);
                        eval(['saveas(h,''EV_THETA_ICA_',filenam,'_dec',int2str(dec),'.fig'')']);
                        ffilenam = [filenam,' EV THETA'];
                        fprf = [ffilenam,' ', num2str(dec)];
                        fprintf(fid,'%s %s\n',fprf);
                        fprintf(fid,'\n',fprf);
                        fprf = [num2str(intraburstivvar(dec)),' ',num2str(extraburstivvar(dec)),' ',num2str(interburstivvar(dec)),' ',...
                                num2str(firstspikevar(dec)),' ',num2str(intraburstivvar_norm(dec)),' ',num2str(extraburstivvar_norm(dec)),...
                                ' ',num2str(interburstivvar_norm(dec)),' ',num2str(firstspikevar_norm(dec))];
                        fprintf(fid,'%s %s %s %s %s %s %s %s\n',fprf);
                        fprintf(fid,'\n',fprf);   
                    end;
                end;
                h2 = figure;
                plot([1:d],intraburstivvar_norm);
                xlim([1 d]);
                x_lim = xlim;
                y_lim = ylim;
                kk = (y_lim(2) - y_lim(1)) / 3;
                if length(vdisc) < 100,
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','r');
                else
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','k');
                end;
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' INTRABURSTIVVAR']);
                    eval(['saveas(h2,''NO_THETA_INTRABURSTIVVAR_',filenam,'.fig'')']);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' INTRABURSTIVVAR']);
                    eval(['saveas(h2,''SP_THETA_INTRABURSTIVVAR_',filenam,'.fig'')']);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' INTRABURSTIVVAR']);
                    eval(['saveas(h2,''EV_THETA_INTRABURSTIVVAR_',filenam,'.fig'')']);
                end;
                h3 = figure;
                plot([1:d],extraburstivvar_norm);
                xlim([1 d]);
                x_lim = xlim;
                y_lim = ylim;
                kk = (y_lim(2) - y_lim(1)) / 3;
                if length(vdisc) < 100,
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','r');
                else
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','k');
                end;
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' EXTRABURSTIVVAR']);
                    eval(['saveas(h3,''NO_THETA_EXTRABURSTIVVAR_',filenam,'.fig'')']);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' EXTRABURSTIVVAR']);
                    eval(['saveas(h3,''SP_THETA_EXTRABURSTIVVAR_',filenam,'.fig'')']);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' EXTRABURSTIVVAR']);
                    eval(['saveas(h3,''EV_THETA_EXTRABURSTIVVAR_',filenam,'.fig'')']);
                end;
                h4 = figure;
                plot([1:d],interburstivvar_norm);
                xlim([1 d]);
                x_lim = xlim;
                y_lim = ylim;
                kk = (y_lim(2) - y_lim(1)) / 3;
                if length(vdisc) < 100,
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','r');
                else
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','k');
                end;
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' INTERBURSTIVVAR']);
                    eval(['saveas(h4,''NO_THETA_INTERBURSTIVVAR_',filenam,'.fig'')']);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' INTERBURSTIVVAR']);
                    eval(['saveas(h4,''SP_THETA_INTERBURSTIVVAR_',filenam,'.fig'')']);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' INTERBURSTIVVAR']);
                    eval(['saveas(h4,''EV_THETA_INTERBURSTIVVAR_',filenam,'.fig'')']);
                end;
                h5 = figure;
                plot([1:d],firstspikevar_norm);
                xlim([1 d]);
                x_lim = xlim;
                y_lim = ylim;
                kk = (y_lim(2) - y_lim(1)) / 3;
                if length(vdisc) < 100,
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','r');
                else
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','k');
                end;
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' FIRSTSPIKEVAR']);
                    eval(['saveas(h5,''NO_THETA_FIRSTSPIKEVAR_',filenam,'.fig'')']);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' FIRSTSPIKEVAR']);
                    eval(['saveas(h5,''SP_THETA_FIRSTSPIKEVAR_',filenam,'.fig'')']);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' FIRSTSPIKEVAR']);
                    eval(['saveas(h5,''EV_THETA_FIRSTSPIKEVAR_',filenam,'.fig'')']);
                end;
                h6 = figure;
                hold on;
                plot([1:d],intraburstivvar_norm,'g');
                plot([1:d],extraburstivvar_norm,'b');
                plot([1:d],interburstivvar_norm,'c');
                plot([1:d],firstspikevar_norm,'m');
                legend('intraburstivvar','extraburstivvar','interburstivvar','firstspikevar',2);
                xlim([1 d]);
                x_lim = xlim;
                y_lim = ylim;
                kk = (y_lim(2) - y_lim(1)) / 3;
                if length(vdisc) < 100,
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','r');
                else
                    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length(vdisc))],'Color','k');
                end;
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),' ALLVAR']);
                    eval(['saveas(h6,''NO_THETA_ALLVAR_',filenam,'.fig'')']);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),' ALLVAR']);
                    eval(['saveas(h6,''SP_THETA_ALLVAR_',filenam,'.fig'')']);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),' ALLVAR']);
                    eval(['saveas(h6,''EV_THETA_ALLVAR_',filenam,'.fig'')']);
                end;
            end;
        end;
    end;
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
fclose(fid);
cd(mmm);