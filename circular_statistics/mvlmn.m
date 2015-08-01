function [ftm, ang, mvl] = mvlmn(ng,degrad)
%MVLMN   Circular mean and mean resultant length.
%   [F, M, R] = MVLMN(NG,DEGRAD) calculates the first trigonometric moment
%   (F), the circular mean (M) and the mean resultant length (R) of a given
%   circular sample (NG). The second input argument (DEGRAD) should contain
%   the information whether the sample is measured in degress ('deg') or
%   radians ('rad'). M is returned in the same measure.
%
%   See also MVL, CIRCULAR_MEAN3 and CIRCULAR_SE.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argument check
error(nargchk(2,2,nargin))

% Main
if isequal(degrad,'deg')    % convert to radians
    ng = ng / 180 * pi;
end
ftm = sum(exp(1).^(1i*ng)) / length(ng);    % first trigonometric moment
ang = angle(ftm);   % mean angle
mvl = abs(ftm);     % mean resultant length
if isequal(degrad,'deg')    % convert to degrees
    ang = ang * 180 / pi;
end