%STORE  Contains useful parts for MATLAB programming.

%Searching for the matching datinx file
    cl = fname(1:6);
    k = dir(pathname);
    j = zeros(length(k),1);
    cell_j = cell(1,length(k));
    for m = 3:length(k),
        cell_j{m} = k(m).name(1:6);
    end;
    j = strcmp(cell_j,cl);
    if isempty(find(j)) == 0,
        filename = k(find(j)).name;
        pont = findstr(filename,'.');
        filenam = filename(1:pont(1)-1);
    end;
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
    
% -------------------------------------------------------------------------------------
% A part of ICA_WAVE
% Wave segmentation regarding the burst localizations
            fdi = find(diff(intraburstivvar_norm));
            lenfdi = length(fdi);
            vbfew = cell(1,lenfdi);
            flag = cell(1,lenfdi);
            burstwave = cell(1,lenfdi);
            linbw = cell(1,lenfdi);
            frames = cell(1,lenfdi);
            contobj = cell(1,lenfdi);
            for b = 2:lenfdi    %first two elements of intraburstivvar_norm are zeros
                vbfew{b} = vdisc(Burst{fdi(b)});    %vbfew contains the loc. of the first skipes of bursts in its first line
                                                    %and loc. of the last spikes of bursts in its second line considering
                                                    %only those steps where Burst matrix changes
                burstlength = vbfew{b}(2,:) - vbfew{b}(1,:) + 1;
                flag{b} = min(burstlength,1000);
                lenvbfew = length(vbfew{2});
                bw = cell(1,lenvbfew);
                lin1 = cell(1,lenvbfew);
                lin2 = cell(1,lenvbfew);
                fr = cell(1,lenvbfew);
                linbw{b} = [];
                frames{b} = [];
                u141 = size(wave,1);
                space = zeros(u141,10) * inf;
                frspace = ones(u141,10);
                for v = 1:lenvbfew
                    lfsp = floor(vbfew{b}(1,v)/newstep);
                    llsp = ceil(vbfew{b}(2,v)/newstep);
                    ff = floor(flag{b}(v)/newstep);                    
                    cf = ceil(flag{b}(v)/newstep);
                    index1 = max(1,lfsp-ff) + (p - 1) * (segmlen / 25);
                    index2 = min(lend,llsp+cf) + (p - 1) * (segmlen / 25);
                    bw{v} = wave(:,index1:index2);    %bw contains a three times longer (but max 2*1000) part of wave for each burst
                                                      %with the burst time in the center
                    if ~isempty(bw{v})
                        linbw{b} = [linbw{b} bw{v} space];    %linbw contains the bw matrices linearly arter each other at every b index
                        lbw = size(bw{v},2);    %creating the frames
                        lf1 = lfsp - index1;
                        lf2 = index2 - llsp;
                        lc = lbw - lf1 - lf2;
                        if (lf1 - 5) > 0 & (lf2 - 5) > 0
                            fr1 = ones(u141,lf1-3);
                            fr21 = ones(5,3);
                            fr22 = zeros(u141-10,3) * inf;
                            fr2 = [fr21; fr22; fr21];
                            fr31 = ones(5,lc);
                            fr32 = zeros(1,lc) * inf;
                            fr33 = ones(u141-12,lc);
                            fr3 = [fr31; fr32; fr33; fr32; fr31];
                            fr4 = fr2;
                            fr5 = ones(u141,lf2-3);
                            fr{v} = [fr1 fr2 fr3 fr4 fr5];
                            frames{b} = [frames{b} fr{v} frspace];
                        else
                            fr{v} = ones(size(bw{v}));
                            frames{b} = [frames{b} fr{v} frspace];
                        end
                        lin1{v} = [2 u141-2 u141-2 2 2];
                        lin2{v} = [lfsp lfsp llsp llsp lfsp];
                    end
                end
                burstwave{b} = bw;  %burstwave contains a cell array at every b index with one matrix per burst cut from wave
                
% Contour plot
%                 levels=2.^[-4:0.1:7];
                contobj{b} = linbw{b} .* frames{b};
                levels = [0:1:10];
                C = contourf_d(abs(contobj{b}),levels);
                figure;contourf_d(abs(linbw{b}),levels);
                figure;contourf_d(abs(linbw{b}),levels);
                for vv = 1:lenvbfew
                    line(lin2{vv},lin1{vv},'Color','k')
                end
                return
            end
            
            
% -------------------------------------------------------------------------------------
% Another part of ICA_WAVE
if p ~= 6
    datinx1 = (p-1) * segmlen + 1; %first point of the interval
    datinx2 = p * segmlen; %last point of the interval
else
    datinx1 = (p-1) * segmlen + 1; %first point of the interval
    datinx2 = lend; %last point of the interval
end

% -------------------------------------------------------------------------------------
% ICA_GUI_2B: Searching for the matching THETA_ICA file
k = dir(pth);
j = zeros(length(k),1);
cell_j = cell(1,length(k));
for m = 3:length(k),
    cell_j{m} = k(m).name(11:30);
end;
j = strcmp(cell_j,cl);
if isempty(find(j)) == 0,
    fln = k(find(j)).name;
end;

% -------------------------------------------------------------------------------------
% ICA_GUI_2B: another way of setting XTickLabel
for fnd = 1:5
    str = ['k_',num2str(fnd),'=f_',num2str(fnd),'*100/pp3;'];
    eval(str);
end
xtl = get(handles.axes1,'XTickLabel');
xtl2 = str2num(xtl);
xtl3 = xtl2 * 100 / pp3;
xtl4 = num2str(xtl3);
set(handles.axes1,'XTickLabel',xtl4)

% -------------------------------------------------------------------------------------
% BACKUP: "last week"
current_date = date;
last_week = cell(1,7);
last_week{1} = current_date;
for c = 1:6
    last_week{c+1} = datestr(datenum(current_date)-c);
end

% -------------------------------------------------------------------------------------
% A part of the thetaselectors: makes disjuct intervals from a set of points:
% points following each other continuusly form an interval, ie. all data points
% of the interval are in the set, and non of the data points out of the intervals
% are in the set
function draw(ip)
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
prepa = zeros(2,lenfdr+1);
prepa(1,1) = ip(1);
for t = 1:lenfdr
    prepa(2,t) = ip(fdr(t));
    prepa(1,t+1) = ip(fdr(t)+1);
end
prepa(2,end) = ip(end);
zs = zeros(1,lenfdr+1);
os = ones(1,lenfdr+1);
patcha = [prepa(1,:); prepa(1,:); prepa(2,:); prepa(2,:); prepa(1,:)];
patchb = [zs; os; os; zs; zs];
ptch = patch(patcha,patchb,'blue');
set(ptch,{'edgecolor'},get(ptch,{'facecolor'}));

% --------------------------------------------------------------------
% From EMPTY_GUI2004: imports base workspace variables
function varargout = edit1_Callback(h, eventdata, handles, varargin)
handles.InputArg = get(handles.edit1,'String');
fstr = findstr(handles.InputArg,',');
fstr(end+1) = length(handles.InputArg) + 1;
lfstr = length(fstr);
ias = cell(1,lfstr);
ia = cell(1,lfstr);
next = 1;
for i = 1:lfstr
    ias{i} = handles.InputArg(next:fstr(i)-1);
    next = fstr(i) + 1;
    cmnd = ['assignin(''caller'',''pia'',' ias{i} ');'];
    evalin('base',cmnd)
    ia{i} = pia;
end
handles.InputArg = ia;
guidata(h,handles);

% -------------------------------------------------------------------------
% From AZSHIFTRUN5: ZMaxLoc histogram.
H3 = figure;
[ZHist,ZHx] = hist(ZMaxLoc,length(ZMaxLoc)/10);
bar(ZHx,ZHist)
cd jpg
eval(['saveas(H3,''',filenam,'_ZHIST.fig'')']);
cd ..
close(H3)
cd mat
str = [filenam '_ZHIST'];
save(str,'ZHist','ZHx')
cd ..