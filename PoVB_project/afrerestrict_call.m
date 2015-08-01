function afrerestrict_call
%AFRERESTRICT_CALL    Calls AFRE_RESTRICT for a sequence of directories.
%
%   See also AFRE_RESTRICT.

% Directories
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl_control\';
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
    afre_restrict(inpdir)
end