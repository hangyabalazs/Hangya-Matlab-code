function czplaceanalysis2_call
%CZPLACEANALYSIS2_CALL    Calls CZPLACEANALYSIS2 for a sequence of directories.
%
%   See also CZPLACEANALYSIS2.

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
for k = 1:length(inpadd)
    fname = inpadd{k};
    czplaceanalysis2b(fname)
end