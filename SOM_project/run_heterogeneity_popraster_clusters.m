load('clustered_cellids')
NumClust = size(clustered_cellids,2);

basedir = ['/Users/ranades/Documents/MATLAB/mPFC paper figures/figures/fig3'];
printfilename = [basedir filesep 'CCGResponse_clusters' datestr(today,'yyyymmdd') '.eps'];

TOTALRESPMATRIX = [];
for iC = 1:NumClust,
    
    % Get cross correlation modulation.
    RESPMATRIX_z = fig3_heterogeneity_popraster(clustered_cellids{iC});
    title(iC)
    if iC == 1,
        MEANRESPMATRIX = nan(NumClust,size(RESPMATRIX_z,2));
    end
    MEANRESPMATRIX(iC,:) = mean(RESPMATRIX_z,1);
    TOTALRESPMATRIX = [TOTALRESPMATRIX ;RESPMATRIX_z];
    % Save figure.
    if iC == 1,
        print('-dpsc2','-r100',printfilename)
    else
        print('-dpsc2','-append','-r100',printfilename)
    end
end

figure
imagesc(MEANRESPMATRIX);
        