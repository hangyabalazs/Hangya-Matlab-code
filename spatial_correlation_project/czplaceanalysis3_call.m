function czplaceanalysis3_call
%CZPLACEANALYSIS3_CALL    Calls CZPLACEANALYSIS3 for a sequence of directories.
%
%   See also CZPLACEANALYSIS3.

% Directories
global DATAPATH
inproot = [DATAPATH 'Czurko\discriminated2\neg\'];
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call APHASE_CORR2B
for k = 1:1
    fname = inpadd{k}
    czplaceanalysis3b(fname)
end