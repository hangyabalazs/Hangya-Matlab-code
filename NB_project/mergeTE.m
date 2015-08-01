fnm = fieldnames(TEa);
TE = TEa;
for k = 1:length(fnm)
    TE.(fnm{k}) = [TE.(fnm{k}) TEb.(fnm{k})];
end

save([fullpth filesep 'TE.mat'],'-struct','TE')