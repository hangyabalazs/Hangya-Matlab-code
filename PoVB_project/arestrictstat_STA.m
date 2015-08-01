function arestrictstat_STA
%ARESTRICTSTAT_STA   Burst statistics for bounded frequency band.
%   ARESTRICTSTAT_STA calculates summary Spike Triggered Average statistics
%   for Po, VPM, PoVPM, LD and nRT data. It plots and saves result figures:
%   one displaying STA indices for all cell and another representing STA
%   indeces as boxplots. Edit code to modify input and output directories!
%
%   See also ARESTRICTSTAT_BURST and AFRE_RESTRICT.

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_stand\'];   % burst analysis data
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_stand\'];
mm = pwd;
dr = dir(inpdir);

% Main
si1_as = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si1_fs = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si1_sp = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si1_afsp = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si2_as = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si2_fs = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si2_sp = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
si2_afsp = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
for o = 3:length(dr)
    inpadd = dr(o).name;    % load phase data
    cd(inpdir)
    cd(inpadd)
    cd('bas')
    ddr = dir(pwd);
    fn = ddr(end).name(1:end-4);
    cmps = strread(fn,'%s','delimiter','_');
    fname = [cmps{1} '_' cmps{2}];
    ff = [inpdir2 fn '_STA.mat'];
    try
        load(ff)
    catch
        lasterr
        continue
    end
    ff = [tabledir 'tablazat_Balazsnak'];   % load position data
    [tbl0 tbl] = xlsread(ff);
    inx = find(strcmp({tbl{:,1}},fname));
    loc = tbl{inx,3};
    
    eval(['si1_as.' loc '(end+1) = sta_index1_as;']);  % STA indices
    eval(['si1_fs.' loc '(end+1) = sta_index1_fs;']);
    eval(['si1_sp.' loc '(end+1) = sta_index1_sp;']);
    eval(['si1_afsp.' loc '(end+1) = sta_index1_afsp;']);
    eval(['si2_as.' loc '(end+1) = sta_index2_as;']);
    eval(['si2_fs.' loc '(end+1) = sta_index2_fs;']);
    eval(['si2_sp.' loc '(end+1) = sta_index2_sp;']);
    eval(['si2_afsp.' loc '(end+1) = sta_index2_afsp;']);
end

% Plot and save
dbclear if error
cd(resdir)

H = figure;
hold on
si1 = si1_as;
clr = 'blue';
pcl = 'b.';
inc = 0;
plot(1+inc,mean(si1.Po),pcl)
plot(2+inc,mean(si1.VB),pcl)
plot(3+inc,mean(si1.PoVPM),pcl)
plot(4+inc,mean(si1.LD),pcl)
plot(5+inc,mean(si1.nRT),pcl)
L1 = line([1+inc 1+inc],[mean(si1.Po)-std(si1.Po) mean(si1.Po)+std(si1.Po)],'Color',clr);
line([2+inc 2+inc],[mean(si1.VB)-std(si1.VB) mean(si1.VB)+std(si1.VB)],'Color',clr')
line([3+inc 3+inc],[mean(si1.PoVPM)-std(si1.PoVPM) mean(si1.PoVPM)+std(si1.PoVPM)],'Color',clr)
line([4+inc 4+inc],[mean(si1.LD)-std(si1.LD) mean(si1.LD)+std(si1.LD)],'Color',clr)
line([5+inc 5+inc],[mean(si1.nRT)-std(si1.nRT) mean(si1.nRT)+std(si1.nRT)],'Color',clr)
si1 = si1_fs;
clr = 'green';
pcl = 'g.';
inc = 0.1;
plot(1+inc,mean(si1.Po),pcl)
plot(2+inc,mean(si1.VB),pcl)
plot(3+inc,mean(si1.PoVPM),pcl)
plot(4+inc,mean(si1.LD),pcl)
plot(5+inc,mean(si1.nRT),pcl)
L2 = line([1+inc 1+inc],[mean(si1.Po)-std(si1.Po) mean(si1.Po)+std(si1.Po)],'Color',clr);
line([2+inc 2+inc],[mean(si1.VB)-std(si1.VB) mean(si1.VB)+std(si1.VB)],'Color',clr')
line([3+inc 3+inc],[mean(si1.PoVPM)-std(si1.PoVPM) mean(si1.PoVPM)+std(si1.PoVPM)],'Color',clr)
line([4+inc 4+inc],[mean(si1.LD)-std(si1.LD) mean(si1.LD)+std(si1.LD)],'Color',clr)
line([5+inc 5+inc],[mean(si1.nRT)-std(si1.nRT) mean(si1.nRT)+std(si1.nRT)],'Color',clr)
si1 = si1_sp;
clr = 'black';
pcl = 'k.';
inc = 0.2;
plot(1+inc,mean(si1.Po),pcl)
plot(2+inc,mean(si1.VB),pcl)
plot(3+inc,mean(si1.PoVPM),pcl)
plot(4+inc,mean(si1.LD),pcl)
plot(5+inc,mean(si1.nRT),pcl)
L3 = line([1+inc 1+inc],[mean(si1.Po)-std(si1.Po) mean(si1.Po)+std(si1.Po)],'Color',clr);
line([2+inc 2+inc],[mean(si1.VB)-std(si1.VB) mean(si1.VB)+std(si1.VB)],'Color',clr')
line([3+inc 3+inc],[mean(si1.PoVPM)-std(si1.PoVPM) mean(si1.PoVPM)+std(si1.PoVPM)],'Color',clr)
line([4+inc 4+inc],[mean(si1.LD)-std(si1.LD) mean(si1.LD)+std(si1.LD)],'Color',clr)
line([5+inc 5+inc],[mean(si1.nRT)-std(si1.nRT) mean(si1.nRT)+std(si1.nRT)],'Color',clr)
si1 = si1_afsp;
clr = 'red';
pcl = 'r.';
inc = 0.3;
plot(1+inc,mean(si1.Po),pcl)
plot(2+inc,mean(si1.VB),pcl)
plot(3+inc,mean(si1.PoVPM),pcl)
plot(4+inc,mean(si1.LD),pcl)
plot(5+inc,mean(si1.nRT),pcl)
L4 = line([1+inc 1+inc],[mean(si1.Po)-std(si1.Po) mean(si1.Po)+std(si1.Po)],'Color',clr);
line([2+inc 2+inc],[mean(si1.VB)-std(si1.VB) mean(si1.VB)+std(si1.VB)],'Color',clr');
line([3+inc 3+inc],[mean(si1.PoVPM)-std(si1.PoVPM) mean(si1.PoVPM)+std(si1.PoVPM)],'Color',clr);
line([4+inc 4+inc],[mean(si1.LD)-std(si1.LD) mean(si1.LD)+std(si1.LD)],'Color',clr);
line([5+inc 5+inc],[mean(si1.nRT)-std(si1.nRT) mean(si1.nRT)+std(si1.nRT)],'Color',clr);
xlim([0 6])
set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'Po','VB','PoVPM','LD','nRT'})
legend([L1 L2 L3 L4],{'all spikes','burst first spikes','single spikes','all first spikes'})
saveas(H,'STA.fig');
saveas(H,'STA.eps');


H = stabox(si1_afsp);
title('STA all first spike')
saveas(H,'STAbox_afsp.fig');
saveas(H,'STAbox_afsp.eps');

cd(mm)

% -------------------------------------------------------------------------
function H = stabox(B)

m1 = B.Po;
m2 = B.VB;
m3 = B.PoVPM;
m4 = B.LD;
m5 = B.nRT;

H = figure;
boxplot([m1 m2 m3 m4 m5],[zeros(size(m1)) ones(size(m2)) 2*ones(size(m3))...
    3*ones(size(m4)) 4*ones(size(m5))],'labels',[{'Po'} {'VB'} {'PoVPM'}...
    {'LD'} {'nRT'}]);
[Wp_eu,Wh_eu] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')