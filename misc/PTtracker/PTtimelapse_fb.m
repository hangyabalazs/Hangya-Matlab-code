function PTtimelapse_fb

% Directoris
imagedir = 'c:\Balazs\_personal\PTtracker\images\';   % images
datadir = 'c:\Balazs\_personal\PTtracker\data\';   % angle data

% File names
pdr = dir(imagedir);
dr = {pdr(3:end).name};
pord = {pdr(3:end).date};  % date of image
pord2 = cellfun(@datenum,pord);
[~, ord] = sort(pord2);   % order by date

% Make frames
next = 1;
NumFrames = length(dr) * 3;
[extension flexion dnm] = deal(nan(NumFrames,1));
for k = 1:2:length(dr)
    inxs = next:next+5;  % frame indices
    open([imagedir dr{ord(k)}])   % open image - extension
    A1 = gca;
    ln = findobj(allchild(A1),'Type','line');   % register to elbow
    xd = get(ln,'XData');
    yd = get(ln,'YData');
    fixpoint = [yd(2) xd(2)];
    I = findobj(allchild(A1),'Type','image');
    cdt = get(I,'CData');
    cdtt = graypad(cdt,300);
    cdtn = cdtt(fixpoint(1):fixpoint(1)+600,fixpoint(2):fixpoint(2)+600,:);
    T = findobj(allchild(A1),'Type','text','Color',[0.75 0 0]);
%     postxt = get(T,'Position');
    set(I,'CData',cdtn);
    set(ln,'YData',yd-fixpoint(1)+300)
    set(ln,'XData',xd-fixpoint(2)+300)
    set(T,'Position',[350 300 0])
    
    open([imagedir dr{ord(k+1)}])   % open image - flexion
    A2 = gca;
    ln = findobj(allchild(A2),'Type','line');   % register to elbow
    xd = get(ln,'XData');
    yd = get(ln,'YData');
    fixpoint = [yd(2) xd(2)];
    I = findobj(allchild(A2),'Type','image');
    cdt = get(I,'CData');
    cdtt = graypad(cdt,300);
    cdtn = cdtt(fixpoint(1):fixpoint(1)+600,fixpoint(2):fixpoint(2)+600,:);
    T = findobj(allchild(A2),'Type','text','Color',[0.75 0 0]);
%     postxt = get(T,'Position');
    set(I,'CData',cdtn);
    set(ln,'YData',yd-fixpoint(1)+300)
    set(ln,'XData',xd-fixpoint(2)+300)
    set(T,'Position',[350 300 0])
    
    H = figure('Position',[555 180 1171 740]);
    A1n = copyobj(A1,H);   % copy extension
    set(A1n,'Units','pixels')
    set(A1n,'Position',[-151.2300 59.7664 887.9693 414.9432])
    
    A2n = copyobj(A2,H);   % copy flexion
    set(A2n,'Units','pixels')
    set(A2n,'Position',[434.2700 59.7664 887.9693 414.9432])
    
    dnm(inxs) = pord2(ord(k));
    text(500,500,datestr(dnm(inxs(end)),1))   % date
    
    A3n = axes('Unit','pixels','Position',[45 510 490 190]);
    savedate = dr{ord(k)}(4:end-4);
    adt = load([datadir 'extension' savedate '.mat']);
    extension(inxs) = adt.angle;
    plot(dnm-pord2(ord(1)),extension,'LineWidth',3)
    ylim([130 175])
    xlim([0 pord2(ord(end))-pord2(ord(1))])
%     set(A3n,'XLimMode','manual','YLimMode','manual')
    line(xlim,[170 170],'Color','k','LineStyle',':')
    axis off
    
    A4n = axes('Unit','pixels','Position',[608 510 490 190]);
    savedate = dr{ord(k+1)}(4:end-4);
    adt = load([datadir 'flexion' savedate '.mat']);
    flexion(inxs) = adt.angle;
    plot(dnm-pord2(ord(1)),flexion,'LineWidth',3)
    ylim([25 50])
    xlim([0 pord2(ord(end))-pord2(ord(1))])
%     set(A4n,'XLimMode','manual','YLimMode','manual')
    line(xlim,[30 30],'Color','k','LineStyle',':')
    axis off
    
    M(inxs) = getframe(H);   % make frame
    close all
    next = next + 6;  % counter
end
M(next:next+100) = M(end);   % freeze last frame

% Make movie
movie2avi(M,'timelapse20.avi','Compression','none','fps',24)

% -------------------------------------------------------------------------
function Xnp = graypad(X,rd)

% Padding
[s1 s2 s3] = size(X);
nans1 = repmat(ones(2*rd+s1,rd)*204,[1 1 s3]);
nans2 = repmat(ones(rd,s2)*204,[1 1 s3]);
Xnp = [nans1 [nans2; X; nans2] nans1];