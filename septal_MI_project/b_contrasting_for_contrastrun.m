function [hh] = b_contrasting_for_contrastrun(vdisc,unit,delta)
%CONTRASTING_FOR_CONTRASTRUN Version of CONTRASTING used by CONTRASTRUN.
%   The functions output argument contains the 'contrasting figure' handle(s).
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