function b_filtertest(c1,c2,s)
%FILTERTEST Tests the effect of filtering on the minimum locations.
%   FILTERTEST(C1,C2,S) refilters the sample S times, from the cutoff frequency
%   C1 to the cutoff frequency C2. C1 and C2 are numbers between 0 and 1, where
%   1 corresponds to the Nyquist frequency. The original filtered data is green
%   on the plots (with cutoff frequency at 7 Hz). The second plot shows only the
%   minimum locations.
%
%   FILTERTEST, by itself, uses C1 = 0.0007 (3.5 Hz), C2 = 0.0035 (17.5 Hz) and
%   S = 5 default values.
%
%   FILTERTEST(S) refilters the sample S times and uses C1 = 0.0007 (3.5 Hz) and 
%   C2 = 0.0035 (17.5 Hz) default values.
%
%   FILTERTEST(C1,C2) refilters the sample 5 times from the cutoff frequency C1
%   to the cutoff frequency C2.
%
%   This function uses FILTFILT zero-phase forward and reverse digital filtering
%   function with FIR1 lowpass filter.
%
%   See also FIR1 and FILTFILT.

% Input arguments check
error(nargchk(0,3,nargin)); 
switch nargin
case 0
    c1 = 0.0007;
    s = 5;
    c2 = 0.0035;
case 1
    s = c1;
    c1 = 0.0007;
    c2 = 0.0035;
case 2
    s = 5;
end;

% Wavephase
[t3,theta,eegs,time] = b_wavephase_for_filtertest;
if isempty(t3),
    return;
end;
figure;

% Finds points halfway between the minimum locations
t4 = diff(t3);
t5 = t4 ./ 2;
t6 = t3(1:(length(t3)-1));
t7 = t6 + t5;
t8 = zeros(1,(length(t7)+2));
t8(1) = 1;
t8(2:(length(t7)+1)) = t7(1:(length(t7)));
t8(length(t7)+2) = length(theta);
t9 = round(t8);

% Refiltering
if length(eegs) > 6150,
    fo = 2048;
else fo = fix(length(eegs)/3);
end;
step = (c2 - c1) / (s - 1);
theta_cell = cell(1,s);
z_cell = cell(1,s);
difs_cell = cell(1,s);
next = 1;
for j = c1:step:c2,
    b = fir1(fo,j);
    theta_cell{next} = filtfilt(b,1,eegs);
    
% Localization of new minimums
    u = zeros(1,(length(t9)-1));
    v = zeros(1,(length(t9)-1));
    z_cell{next} = zeros(1,(length(t9)-1));
    for i = 1:(length(t9)-1),
        u(i) = min(theta_cell{next}(t9(i):t9(i+1)));
        v(i) = find(theta_cell{next}(t9(i):t9(i+1))==u(i));
        z_cell{next}(i) = t9(i) + v(i) - 1;
    end;

% Difference between the minimum locations
    difs_cell{next} = zeros(1,(length(t9)-1));
    difs_cell{next} = abs(t3-z_cell{next});
    plot(theta_cell{next});
    hold on
    next = next + 1;
end;
plot(theta,'g');

% Plotting lines at the minimum localizations
figure;
for i = 1:s,
    h = length(z_cell{i});
    for j = 1:h,
    rajz = [z_cell{i}(j) z_cell{i}(j);-2 2];
    line(rajz(1,:),rajz(2,:));
    text(rajz(1,1),rajz(2,1)+1+i/5,int2str(i),'FontSize',8,...
            'HorizontalAlignment','center');
    end;
end;
for j = 1:length(t3),
    rajz = [time(t3(j)) time(t3(j));-2 2];
    L = line(rajz(1,:),rajz(2,:));
    set(L,'Color','g');
end;