function vs = createvs(v,la)
%CREATEVS   Creates supressed data matrices for ANALYSE GUI.
%   ANALYSE GUI, in order to speed up the display of EEG and unit, saves
%   supressed data files ('vs') at first plotting. During oncoming openings,
%   it loads first only 'vs' matrices. VS = CREATEVS(V,LA) needs EEG or unit
%   in V, the base of logarithm used for segmenitng V (LA), and returns the
%   supressed data in VS.
%
%   See also ANALYSE and LINER.

i = 1;
vv = v;
vs = cell(1,floor(log(length(v))/log(la)));
while length(vv) >= 2*la   
    vs{i} = vv;
    vv = [vv,zeros(1,2*la*ceil(length(vv)/2/la)-length(vv))];
    vv = reshape(vv,2*la,length(vv)/2/la);
    minvv = min(vv,[],1); 
    maxvv = max(vv,[],1);
    vv = reshape([minvv;maxvv],1,length(minvv)*2);
    i = i+1;
end
vs{i} = vv;