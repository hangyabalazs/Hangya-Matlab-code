function list = b_natclusttest
%NATCLUSTTEST   Search for reasonable cutting ICC value.
%   NATCLUST runs on a sequence of files and determines the range of cutting ICC
%   (inconsistency coefficient) value, with whom the number of clusters by natural
%   cluster analysing strategy falls between 2 and 20. You can modify the input
%   directory through editing the program code.
%   The result is given in a structure array with fields 'name' - filename, 'limits'
%   - segment boundaries, 'method' - clustering method ('average', 'single' or 'ward')
%   and 'range' - the ICC values. The result is saved in the output directory which
%   you can also  modify in the progrtam code.
%
%   See also NATCLUST and NATCLUSTTEST2.

% Input argument check
error(nargchk(0,0,nargin));

% Input directory
global DATAPATH
where = [DATAPATH,'DATA\analysenow3\'];    %Here are the data files
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'ICA\natclust\']);  %Here are the results

% Progress indicator
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
list = [];      %initialize list of analysed files
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
        
% Cluster analysis
        for p = 1:pend      %"THETA SEGMENT CYCLE"
            i_frst = out_int(2*p-1);    %create indeces
            i_scnd = out_int(2*p);
            ind1 = i_frst + (datinx1 - 1);
            ind2 = i_scnd + (datinx1 - 1);
            lgt = i_scnd - i_frst;
            
            for lnk = 1:3       %"METHOD CYCLE"
                switch lnk
                case 1
                    method = 'ward';
                case 2
                    method = 'single';
                case 3
                    method = 'average';
                end
                
                % DIST, LINKAGE and COPHENET
                ivs = diff(vd{p});
                livs = length(ivs);
                dmtx = zeros(livs,2);
                iivs = ivs';
                dmtx(:,1) = iivs;
                dist = pdist2(dmtx);
                links = linkage(dist,method);
                
                format long
                    
                ICC = 1.1547;
                clusno = inf;
                while clusno > 20
                    % CLUSTER
                    c = cluster(links,ICC);
                    clusno = max(c)    %number of clusters
                    ICC = ICC + 0.00001;
                end
                I1 = ICC - 0.00001;
                
                ICC = 1.1548;
                clusno = -inf;
                while clusno < 2
                    % CLUSTER
                    c = cluster(links,ICC);
                    clusno = max(c)    %number of clusters
                    ICC = ICC - 0.00001;
                end
                I2 = ICC + 0.00001;
                
                if islist    %generate list of analysed files
                    lenlist = length(list);
                    list(lenlist+1).name = filenam;
                    list(lenlist+1).limits = [ind1 ind2];
                    list(lenlist+1).method = method;
                    list(lenlist+1).range = [I2 I1];
                end
            end     %end of "method cycle"
        end     %end of "theta segment cycle"
    end     %end of "registration segment cycle"
    waitbar(o/sf)   %Progress indicator
end     %end of "cell cycle"
save list list

% Close
close(wb);   %Close progress indicator
close all