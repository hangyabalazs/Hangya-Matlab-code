function ezshift_xlsprocess3_zshift
%EZSHIFT_XLSPROCESS3_ZSHIFT    Generates plots from ZSHIFT output excel file.
%   EZSHIFT_XLSPROCESS3_ZSHIFT works on the output excel file of EZSHIFTRUN.
%   It compares shift values between cells with different MI values.
%
%   See also EZSHIFTRUN, EZSHIFT_XLSPROCESS3 and ENTRYRUN3C_ZSHIFT.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
dr = [DATAPATH 'Ezshift\'];
fn = [dr 'zshift_theta.xls'];
inpdir_theta = [DATAPATH 'MI\Theta\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta\control\'];
mm = pwd;
cd(dr)

% Files
[files_theta, files_short_theta] = genfiles(inpdir_theta);

% Read excel file
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'Sheet1');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Plot I
MIp = [];
MIn = [];
MIna = [];
H1 = figure;
hold on
for i = 1:size(atx,1)
    cd(inpdir_theta)    % load MI values
    inx = find(strcmp(atx{i,1},files_short_theta));
    if isempty(inx)
        Wh_theta = [];
    else
        ftheta = files_theta{inx};
        load(ftheta)
        MI_real_theta = aIxy;
        cd(inpdir_theta_control)
        ftheta = files_theta{inx};
        load(ftheta)
        MI_control_theta = aIxy;
        if median(MI_real_theta) > median(MI_control_theta)     % check if MI > control
            [Wp_theta,Wh_theta] = b_signrank2(MI_real_theta,MI_control_theta,'alpha',0.01);
        else
            Wh_theta = 0;
        end
    end
    if ~isempty(Wh_theta)
        if Wh_theta
            h1 = plot(atx{i,4},atx{i,5},'r.','MarkerSize',20);
            MIp(end+1) = atx{i,4};
        else
            h2 = plot(atx{i,4},atx{i,5},'g.','MarkerSize',20);
            MIn(end+1) = atx{i,4};
        end
    else
        h2 = plot(atx{i,4},atx{i,5},'k.','MarkerSize',20);
        MIna(end+1) = atx{i,4};
    end
end
saveas(H1,'zshift_all_MI.fig')

H2 = mibox(MIp(~isnan(MIp)),MIn(~isnan(MIn)),'MI sig.','MI non-sig.');
saveas(H2,'zshift_all_MIbox.fig')



% -------------------------------------------------------------------------
function H = mibox(m1,m2,l1,l2)

H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{l1} {l2}]);
[Wp,Wh] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')

% -------------------------------------------------------------------------
function [files, files_short] = genfiles(inpdir)

mm = pwd;
cd(inpdir)
files0 = dir(pwd);
files = {};
files_short = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if isequal(ext,'.mat')
        files{end+1} = files0(i).name;
        files_short{end+1} = files0(i).name(1:6);
    end
end
cd(mm)