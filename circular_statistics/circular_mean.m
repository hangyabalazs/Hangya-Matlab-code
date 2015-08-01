function ang = circular_mean(ng,degrad)
%CIRCULAR_MEAN   Circular mean.
%   M = CIRCULAR_MEAN(NG,DEGRAD) calculates the circular mean (M) of a
%   given circular sample (NG). The second input argument (DEGRAD) should
%   contain the information whether the sample is measured in degress
%   ('deg') or radians ('rad'). M is returned in the same measure.
%
%   See also MVL, MVLMN and CIRCULAR_SE.

% Input argument check
error(nargchk(2,2,nargin))

% Main
if isequal(degrad,'deg')    % convert to radian
    ng = ng / 180 * pi;
end
ftm = sum(exp(1).^(i*ng)) / length(ng);    % first trigonometric moment
ang = angle(ftm);   % mean angle
if isequal(degrad,'deg')    % convert to degrees
    ang = ang * 180 / pi;
end