function pnophase
%PNOPHASE   Phase histogram.
%   PNOPHASE calculates and plots phase histogram for spikes of PNO neurons
%   with respect to cortical EEG or multiunit. 'Pizzaplot ' is also plotted
%   between circular mean +- 1 circular SD.
%
%   See also PIZZAPLOT.

% Directories
global DATADIR
global DATAPATH
inpdir1 = [DATADIR 'PViktor\raw\'];
inpdir2 = [DATADIR 'Pviktor\disc\'];
resdir = [DATAPATH 'PViktor\phase2\'];

% Filelist
files = dir(inpdir1);
files = files(3:end);
sf = length(files);

% Main
H1 = figure;
for o = 2:2
    
    % Load raw data
    fname = files(o).name;
    ff = [inpdir1 fname];
    rawdata = load(ff);
    
    % Load discriminated data
    ff2 = [inpdir2 fname(1:end-4) '_d.mat'];
    load(ff2);
    
    % Check EEG sampling rate
    eeg = rawdata.eeg.values';
    eeg = double(eeg);
    sr_eeg = round(1/(rawdata.eeg.times(2)-rawdata.eeg.times(1)));
    if ~isequal(sr_eeg,1000)
        error('EEG should be sampled on 1000 Hz.')
    end
    
    % Check MUA envelope sampling rate
    mua = rawdata.muacover.values';
    mua = double(mua);
    sr_mua = round(1/(rawdata.muacover.times(2)-rawdata.muacover.times(1)));
    if ~isequal(sr_mua,1000)
        error('MUA envelope should be sampled on 1000 Hz.')
    end
    
    % Adjust discriminated unit to EEG sampling rate
    sr_unit = round(1/(rawdata.unit.times(2)-rawdata.unit.times(1)));
    cst = sr_unit / sr_eeg;
    vdisc = round(vdisc/cst);
    sr = 1000;
    
    % Use indicator variable to blank out spikes
    if isfield(rawdata,'level')
        disp(['Levels detected for ' fname '.'])
        tms = [rawdata.level.times; rawdata.eeg.times(end)];   % step times
        lvs = rawdata.level.level;   % levels
        lvind = zeros(1,length(eeg));   % make indicator variable from step times and levels
        ons = find(lvs);
        for k = 1:length(ons)
            inx1 = round(tms(ons(k))*sr) + 1;
            inx2 = round(tms(ons(k)+1)*sr);
            lvind(inx1:inx2) = 1;
        end
        vdisc(~logical(lvind(vdisc))) = [];
    else
        disp(['No levels detected for ' fname '.'])
    end
    
    % Phase calculation - EEG
    [ang, mvl, ~, aang, ~, ~, Hp] = phasehist(eeg,vdisc,sr,0,4);    % PHASE
    close(Hp)
    n = length(aang);
    astd = b_circular_std(aang);
    aang = aang * 180 / pi;
    ang = ang * 180 / pi;
    astd = astd * 180 / pi;
    edges = -180:20:180;     % edges for phase histogram
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nm,xout] = histc(aang,edges);   % phase histogram
    nm = nm(1:end-1);
        
    % Plot
    figure(H1);
    subplot(2,2,1)
    B = bar(cnts,nm'/n);
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    ach = allchild(H1);     % figure title
    ax = findobj(ach,'type','axes');
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
    title(ax(end),titlestr)
    x_lim = xlim;
    y_lim = ylim;
    str = ['\it{Mean angle: }' '\bf ' num2str(ang)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','k')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','k')
    str = ['\it{Circulr SD: }' '\bf ' num2str(astd)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','k')
    str = ['\it{n: }' '\bf ' num2str(n)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/4,str,'Color','k')
    
    subplot(2,2,2)
    pizzaplot(deg2rad(ang-astd),deg2rad(ang+astd),[0.16 0.38 0.27])
    
     % Phase calculation - MUA
    [ang, mvl, ~, aang, ~, ~, Hp] = phasehist(mua,vdisc,sr,0,4);    % PHASE
    close(Hp)
    n = length(aang);
    astd = b_circular_std(aang);
    aang = aang * 180 / pi;
    ang = ang * 180 / pi;
    astd = astd * 180 / pi;
    edges = -180:20:180;     % edges for phase histogram
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nm,xout] = histc(aang,edges);   % phase histogram
    nm = nm(1:end-1);
    
    % Plot
    figure(H1);
    subplot(2,2,3)
    B = bar(cnts,nm'/n);
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    ach = allchild(H1);     % figure title
    ax = findobj(ach,'type','axes');
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
    title(ax(end),titlestr)
    x_lim = xlim;
    y_lim = ylim;
    str = ['\it{Mean angle: }' '\bf ' num2str(ang)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','k')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','k')
    str = ['\it{Circulr SD: }' '\bf ' num2str(astd)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','k')
    str = ['\it{n: }' '\bf ' num2str(n)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/4,str,'Color','k')
    
    subplot(2,2,4)
    pizzaplot(deg2rad(ang-astd),deg2rad(ang+astd),[0.16 0.38 0.27])
    
    maximize_figure(H1)
    
% Save
    fn11 = [resdir fname(1:end-4) '_PHASE.fig'];
    saveas(H1,fn11)
%     fn12 = [resdir fname(1:end-4) '_PHASE.mat'];
%     save(fn12,'aang','ang','mvl')
end
close(H1)