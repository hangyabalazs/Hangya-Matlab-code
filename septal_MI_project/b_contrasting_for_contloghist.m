function [dsc] = b_contrasting_for_contloghist(vdisc,unit,delta)
%CONTRASTING_FOR_CONTLOGHIST Version of CONTRASTING used by CONTRASTLOGHIST.
%
%   The functions output argument contains the contrasted vdisc (discriminated unit).
%
%   See also CONTRASTING, CONTRASTRUN and CONTRASTLOGHIST.

% Input arguments check
error(nargchk(3,3,nargin));

% Contrasting
close all
% figure;
% plot(unit)
% title(['unit']);
% waitforbuttonpress
dsc = zeros(1,length(unit));
dsc(vdisc) = 1;
dsc_rem = dsc;
for i = 1:length(delta),
    if vdisc(2) > vdisc(1) + delta(i)
        dsc(vdisc(1)) = 0;
    end;
    for j = 2:length(vdisc)-1,
        if vdisc(j-1) < vdisc(j) - delta(i) & vdisc(j+1) > vdisc(j) + delta(i),
            dsc(vdisc(j)) = 0;
        end;
    end;
    if vdisc(end-1) < vdisc(end) - delta(i)
        dsc(vdisc(end)) = 0;
    end;
    hh(i) = figure;
    subplot(2,1,1);
    plot(dsc)
    % title(['delta = ',num2str(delta(i))]);
    subplot(2,1,2);
    plot(unit)
    % title(['unit']);
    % waitforbuttonpress
    % dsc = dsc_rem;
end;