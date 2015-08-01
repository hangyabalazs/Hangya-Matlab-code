function b_entryzshift_compare
%ENTRYZSHIFT_COMPARE   Saves combined entropy and zshifted-entropy files.
%   ENTRYZSHIFT_COMPARE saves plots containing uncertainty coefficients and
%   zshifted uncertainty coefficients. It also saves these data in mat
%   files.
%
%   See also ENTRYRUN_ZSHIFT.

% Input argumnet check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir1 = [DATAPATH 'Entry_onset\Theta\entropy\power\line\windowsize1000_overlap1\real\'];
inpdir2 = [DATAPATH 'Entry_zshift_onset\Theta\entropy\power\line\windowsize1000_overlap1\real\'];
resdir = [DATAPATH 'Entry_zshift_onset\Compare\'];
[files files_short files_long] = filelist(inpdir2);
files_short2 = unique(files_short);
sf = length(files_short2);

% Open
for i = 1:sf
    files_short2{i}
    ind = find(strcmp(files_short,files_short2{i}));
    cmps = strread(files(ind(1)).name,'%s','delimiter','_');
    fn = [cmps{1} '_' cmps{2} '_' cmps{3} '_' cmps{4} '_' cmps{5} '_ENTROPY1line.fig'];
    cd(inpdir1)
    open(fn)
    ch1 = allchild(gcf);
    axs1 = findobj(ch1,'Type','axes');
    axes(axs1(find(axs1==min(axs1))))
    chh11 = allchild(gca);
    lns11 = findobj(chh11,'Type','line');
    unit_eeg1 = get(lns11,'YData');
    axes(axs1(find(axs1~=min(axs1)&axs1~=max(axs1))))
    chh12 = allchild(gca);
    lns12 = findobj(chh12,'Type','line');
    eeg_unit1 = get(lns12,'YData');
    
    cd(inpdir2)
    open(fn)
    ch2 = allchild(gcf);
    axs2 = findobj(ch2,'Type','axes');
    axes(axs2(find(axs2==min(axs2))))
    chh21 = allchild(gca);
    lns21 = findobj(chh21,'Type','line');
    unit_eeg2 = get(lns21,'YData');
    axes(axs2(find(axs2~=min(axs2)&axs2~=max(axs2))))
    chh22 = allchild(gca);
    lns22 = findobj(chh22,'Type','line');
    eeg_unit2 = get(lns22,'YData');
    
% Plot
    figure
    hold on
    plot(unit_eeg1,'r')
    plot(eeg_unit1,'b')
    plot(unit_eeg2,'Color',[255 191 0]/256)
    plot(eeg_unit2,'Color',[0 191 255]/256)
    
% Save
    cd(resdir)
    figname = [cmps{1} '_' cmps{2} '_' cmps{3} '_' cmps{4} '_' cmps{5} '_ENTRYZSHIFT.fig'];
    saveas(gcf,figname)
    matname = [cmps{1} '_' cmps{2} '_' cmps{3} '_' cmps{4} '_' cmps{5} '_ENTRYZSHIFT.mat'];
    UnitEeg1 = unit_eeg1;
    UnitEeg2 = unit_eeg2;
    EegUnit1 = eeg_unit1;
    EegUnit2 = eeg_unit2;
    save(matname,'UnitEeg1','UnitEeg2','EegUnit1','EegUnit2')
    
    close all
end    
    

% -------------------------------------------------------------------------
function [files2, files2_short files2_long] = filelist(inpdir)

files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
files2_long = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name(1:min(6,length(files(i).name)));
        files2_long{end+1} = files(i).name;
    end
end
files2 = files2(2:end);