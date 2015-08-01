function b_printdata(list_theta,list_noth)
%PRINTDATA   Prints entropy figures.
%   PRINTDATA(LT,LN) prints raw data and uncertainty coefficients
%   ('ENTROPY3line' fig) for the theta and non-theta segments of the given
%   cells. Input argumnets LT and LN should contain the segments in cells
%   of strings of the following form: 'cell id' '_' 'datinx1 '_' 'datinx2'
%   where 'datinx1' and 'datinx2' are the bounderies of theta or non-theta
%   segment.
%
%   See also ENTRYRUN3.

% Input argumnet check
error(nargchk(2,2,nargin))
if ~isequal(length(list_theta),length(list_noth))
    error('Input argumnets should be of the same length.')
end

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'iontophoresis\'];
entrydir = [DATAPATH 'Entry3\Theta\entropy\power\line\windowsize1000_overlap1\real\'];
resdir = [DATAPATH 'Entry3\Stat\figs\'];
mm = pwd;
cd(inpdir)
[files files_short] = b_filelist(inpdir);

% Import segments
for i = 1:length(list_theta)
    cmps = strread(list_theta{i},'%s','delimiter','_');
    name_theta = [cmps{1} '_' cmps{2}];
    datinx1_theta = str2num(cmps{3});
    datinx2_theta = str2num(cmps{4});
    seglen_theta = datinx2_theta - datinx1_theta;
    cmps = strread(list_noth{i},'%s','delimiter','_');
    name_noth = [cmps{1} '_' cmps{2}];
    datinx1_noth = str2num(cmps{3});
    datinx2_noth = str2num(cmps{4});
    seglen_noth = datinx2_noth - datinx1_noth;
    if ~isequal(name_theta,name_noth)
        error('Technical error 19')
    else
        name = name_theta;
    end
    
    inx = find(strcmp(name,files_short));
    cd(inpdir)
    data = b_load_data(files(inx(1)).name);
    unit = data(:,2);
    unit = unit';
    eeg = data(:,1);
    eeg = eeg';
    unit_theta = unit(datinx1_theta:datinx2_theta);
    unit_noth = unit(datinx1_noth:datinx2_noth);
    unit_new = [unit_theta unit_noth];
    eeg_theta = eeg(datinx1_theta:datinx2_theta);
    eeg_noth = eeg(datinx1_noth:datinx2_noth);
    eeg_new = [eeg_theta eeg_noth];
    seglen_new = seglen_theta + seglen_noth;
    time = linspace(0,seglen_new, length(eeg_new));
    
% Plot
    H = figure;
    S1 = subplot(3,1,1);
    plot(time,eeg_new)
    y_lim = ylim;
    line([seglen_theta seglen_theta],[y_lim(1) y_lim(2)],'Color','green');
    axis([time(1) time(end) y_lim(1) y_lim(2)])
    S2 = subplot(3,1,2);
    plot(time,unit_new)
    y_lim = ylim;
    line([seglen_theta seglen_theta],[y_lim(1) y_lim(2)],'Color','green');
    axis([time(1) time(end) y_lim(1) y_lim(2)])
    S3 = subplot(3,1,3);
    cd(entrydir)
    fnm = [name(1:6) '_ENTROPY3line.fig'];
    E = open(fnm);
    EA = gca;
    AC = allchild(EA);
    bl = findobj(AC,'Color','blue');
    set(bl,'LineStyle','--')
    copyobj(AC,S3)
    close(E)
    axis tight
    subplot(S1)
    title([name(1:3) ' ' name(5:6)])
    
% Save
    cd(resdir)
    saveas(H,name,'fig')
    
% Print
    set(H,'PaperOrientation','landscape','PaperPositionMode','manual',...
        'PaperPosition',[0 0 28 21])
    str = ['print -f' num2str(H)];
    eval(str)
    close all
end