function afrerestrict_gamma_call
%AFRERESTRICT_GAMMA_CALL    Calls AFRE_RESTRICT_GAMMA for a sequence of directories.
%
%   See also AFRE_RESTRICT_GAMMA.

% Directories
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call AFRERESTRICT
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    afre_restrict_gamma(inpdir)
end