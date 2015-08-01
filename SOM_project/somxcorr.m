function somxcorr(SpkData)
%SOMXCORR Cross-correlation for SOM neurons.
%   SOMXCORR(SPKDATA) calls SOMCCG for a sequence of files to calculate
%   original and normalized cross-correlogram (see SOMCCG for details). The
%   input struct SPKDATA should contain two fields for every neruons:
%   'spikes' with action potential timestamps in seconds and 'cellid'.
%
%   See also SOMCCG.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
resdir = [DATAPATH 'SOM\CCG2\fig\'];

% Cross-correlation
cellno = length(SpkData);   % number of cells
wn = 200;    % window size in ms
H = figure;
for k1 = 1:cellno
    for k2 = k1:cellno
        disp([num2str(k1) ' ' num2str(k2)])
        vd1 = SpkData(k1).spikes;    % timestamps for cell1 (in sec.)
        vd2 = SpkData(k2).spikes;    % timestamps for cell2 (in sec.)
        [H1 H2 trsc] = somccg(vd1,vd2,wn);  % CCG
        
        figure(H1)  % modify plots
        A = gca;
        ach = findobj(allchild(A),'Type','hggroup');
        set(ach,'FaceColor',[0 0 0],'BarWidth',1)
        set(A,'XLim',[-30 30])
        axis off
        if ~isempty(H2)
            figure(H2)
            A = gca;
            ach = findobj(allchild(A),'Type','hggroup');
            set(ach,'FaceColor',[0 0 0],'BarWidth',1)
            set(A,'XLim',[-30 30])
            axis off
        end
        
        cell1id = SpkData(k1).cellid{1};   % figure title
        cell1id(cell1id=='.') = '_';
        cell1idt = cell1id;
        cell1idt(cell1idt=='_') = ' ';
        cell2id = SpkData(k2).cellid{1};
        cell2id(cell2id=='.') = '_';
        cell2idt = cell2id;
        cell2idt(cell2idt=='_') = ' ';
        figure(H1)      % summary figure with subplots
        aA1 = allchild(gca);
        figure(H)
        S = subplot(cellno,cellno,(k1-1)*cellno+k2);
        copyobj(aA1,S)
        xlim([-wn wn])
        figure(H1)
        title([cell1idt ' -> ' cell2idt])
        saveas(H1,[resdir 'CCG_' cell1id '_vs_' cell2id '.fig'])    % save
        saveas(H1,[resdir 'CCG_' cell1id '_vs_' cell2id '.tif'])
        saveas(H1,[resdir 'CCG_' cell1id '_vs_' cell2id '.eps'])
        if ~isempty(H2)     % no H2 for auto-correlation
            figure(H2)
            title([cell1idt ' -> ' cell2idt])
            saveas(H2,[resdir 'normCCG_' cell1id '_vs_' cell2id '.fig'])
            saveas(H2,[resdir 'normCCG_' cell1id '_vs_' cell2id '.tif'])
            saveas(H2,[resdir 'normCCG_' cell1id '_vs_' cell2id '.eps'])
        end
        close(H1,H2)
    end
end
hgsave2(H,[resdir 'allCCG_' cell1id '.fig'])    % save
close all