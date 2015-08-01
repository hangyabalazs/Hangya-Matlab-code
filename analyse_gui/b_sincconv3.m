function zint = b_sincconv3(time,unit,newstep,fs)
%SINCCONV3   Convulation with sinc function.
%   Z = SINCCONV(TIME,UNIT,NEWSTEP,FS) convolves the input UNIT with sinc function,
%   when NEWSTEP is the rate between the new and the former sampling frequency (FS)
%   by the resampling of the eeg (it determins the length and the resolution
%   of the output).
%
%   Automatic thresholding is done by an inner thresholding subfunction. Result of
%   thresholding is plotted.
%
%   See also SINC and SINCCONV.

% Thresholding and discrimination
T = newthres(time,unit);
vdisc = disc(unit,T);

% Variable transformation regarding the new sampling frequency
nunit = length(unit);
dt = 1 / fs;
fcut = 100; 
fsnew = fs / newstep;
dtnew = newstep / fs;
fsratio = fsnew / fs;
told = vdisc * dt * fcut;          
tnew = (1:nunit*fsratio) * dtnew * fcut;  

% Convolution
lentold = length(told);
zint = 0;     % final sum
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end

% --------------------------------------------------------------------------
function T = newthres(time,unit)

% Checking if base line is stable
d = 10000;
dd = 1;
mn = [];
while dd < time(end)
    mn(end+1) = mean(unit(dd:dd+d));
    dd = dd + d;
end
mn(end+1) = mean(unit(dd:end));
if std(mn) > 0.1
    str = [fname(1:6) ' : Instable base line.'];
    disp(str)
end

% Action potential amplitudes
unit2 = unit;
for k = 1:10
    m(k) = max(unit2);
    unit2(find(unit2==m(k))) = [];
end
MS = min(m(6:10));

% Mean
MI = mean(unit);

% Threshold
T = MI + 2 * (MS - MI) / 5;

% Plotting
H = figure;
plot(time,unit)
line([time(1) time(end)],[MI MI],'Color','r')
line([time(1) time(end)],[MS MS],'Color','r')
line([time(1) time(end)],[T T],'Color','g')
% scr = get(0,'ScreenSize');
% set(H,'MenuBar','none')
% set(H,'Position',[scr(1)+1 scr(4)-170 300 150])
% drawnow
T2 = [];
T2 = input('Give the threshold! [Press Enter if the orginal threshold is suitable] ');
if ~isempty(T2)
    T = T2;
end