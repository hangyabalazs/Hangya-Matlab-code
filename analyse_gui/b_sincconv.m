function zint = b_sincconv(unit,newstep)
%SINCCONV   Convulation with sinc function.
%   Z = SINCCONV(UNIT,NEWSTEP) convulates the input UNIT with sinc function, when
%   NEWSTEP is the rate between the new and the former sampling frequency, i.e.
%   10 kHz by the resampling of the eeg (it determins the length and the resolution
%   of the output).
%
%   Automatic thresholding is done by an inner thresholding subfunction. Result of
%   thresholding is plotted.
%
%   See also SINC and SINCCONV2.

% Thresholding and discrimination
T = newthres;
b_disc(T)
b_var2ws('disc','caller')

% Variable transformation regarding the new sampling frequency
nunit = length(unit);
fs = 10000;       % old sampling frequency in Hz
dt = 1 / fs;
fcut = 100; 
fsnew = 10000 / newstep;
dtnew = newstep / 10000;
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
function T = newthres

%Reading the input data
global IN
if isempty(IN)
    error('Data has not been inported.')
end
b_var2ws('in','caller')

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
t = ['THRESHOLD ',fname(1:3),' ',fname(5:6),' ',num2str(datinx1),' ',num2str(datinx2)];
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))/5,MS+(y_lim(2)-MS)/2,t,'FontSize',8);
scr = get(0,'ScreenSize');
set(H,'MenuBar','none')
set(H,'Position',[scr(1)+1 scr(4)-170 300 150])
drawnow