function period_handle = b_lombper_for_ica(fname,i_first,i_second,dec)
%LOMBPER_FOR_ICA   Creates Lomb periodogram plots for ICA_GUI2B.
%   LOMBPER_FOR_ICA(FNAME,I_FIRST,I_SECOND,DEC) needs four input arguments: FNAME
%   is the file name, I_FIRST and I_SECOND are the limits of the interval and DEC
%   is the number of clusters (see ICA_WAVE for more information).
%
%   The function loads data previously saved by PERIODRUN, creates the plot and
%   return its handle.
%
%   See also ICA_GUI2B, LOMBPER, LOMB_PERIOD and PERIODRUN.

% Input arguments check
error(nargchk(4,4,nargin));

% Load
global DATAPATH
fnm = ['LOMB_' fname '_' num2str(dec) '_' num2str(i_first) '_' num2str(i_second)];
ff = fullfile(DATAPATH,'ICA\ica_gui2b\lomb',fnm);
load(ff)

z05 = z(1);
z02 = z(2);
z01 = z(3);
z001 = z(4);

% Plot Lomb-spectrum and signif. levels
period_handle = fill(pxx,pyy,'r');
y_lim = ylim;
axis([0 max(pxx) y_lim(1) y_lim(2)])
hold on
v = axis;
% plot([v(1) v(2)], [ z31 z31], 'k:')
% text(v(2),z31,' 31% ')
plot([v(1) v(2)], [ z05 z05], 'k:')
text(v(2),z05,' 5% ')
plot([v(1) v(2)], [ z02 z02], 'k:')
text(v(2)+v(2)*.03,z02,' 2% ')
plot([v(1) v(2)], [ z01 z01], 'k:')
text(v(2)+v(2)*.06,z01,' 1% ')
plot([v(1) v(2)], [ z001 z001], 'k:')
text(v(2),z001,' .1% ')
xlabel(' frequency [Hz]')
ylabel(' PSD [s^2]')
% title(' Lomb - Scargle periodogram ')
hold off