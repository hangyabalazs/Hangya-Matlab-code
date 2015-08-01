function afrerestrictstand_call
%AFRERESTRICTSTAND_CALL    Calls AFRE_RESTRICT_STAND for a sequence of directories.
%
%   See also AFRE_RESTRICT_STAND.

% Directories
inproot = 'Y:\extra32\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 1:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call AFRE_RESTRICT_STAND
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    afre_restrict_stand(inpdir)
end