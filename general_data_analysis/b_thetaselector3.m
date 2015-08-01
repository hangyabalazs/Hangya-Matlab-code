function [H,OM] = b_thetaselector3(power,f,newstep)
%THETASELECTOR3   Theta selection upon wavelet power maximum localization.
%   [H,OM] = THETASELECTOR3(POWER,F,NEWSTEP) find maximum localization of POWER
%   in each coloumn, using F scalevector and NEWSTEP constant (which characterizes
%   downsampling of the EEG sampled on 10 kHz). POWER is plotted on figure H with
%   white pots signaling the maximums. Corrected maximum value (10th high value)
%   is written in the image. A bar shows the theta segments (segments where
%   maximums fall in the theta band).
%
%   OM 2-N output matrix contains the maximum localization for each coloumn (time
%   dimension), and the maximum value itself.
%
%   Theta band: 2.5 - 6 Hz.
%
%   See also THETASELECTOR and THETASELECTORRUN.

% Input arguments check
error(nargchk(3,3,nargin));

% Close all figures
close all

% Finding theta band
fnd = find(f>6);
pwind1 = fnd(end);
fnd = find(f<2.5);
pwind2 = fnd(1);

% Corrected maximum value
sw1 = size(power,1);
sw2 = size(power,2);
time = (([1:sw2] - 1) * newstep + 1) / 10000;
lenw = sw1 * sw2;

power2 = reshape(power,1,lenw);
spower = sort(power2);
mm = spower(end-10);
clear power2 spower

% Maximum localizations
maxes = max(power);
maxloc = zeros(1,sw2);
for i = 1:sw2
    maxloc(i) = find(power(:,i)==maxes(i));
end

% Plotting
H = figure;
subplot('Position',[0.13 0.11+0.1 0.775 0.815-0.1])
imagesc(power)
hold on
plot(maxloc,'Color',[1 1 1],'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 1 1],...
    'MarkerSize',0.5,'LineStyle','none')
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
line([x_lim(1) x_lim(2)],[pwind1 pwind1],'Color','g');
line([x_lim(1) x_lim(2)],[pwind2 pwind2],'Color','g');
xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
text(xcoord,ycoord,num2str(mm),'Color',[1 1 1]);
ff = round(f*100) / 100;
time = round(time*100) / 100;
b_rescalefig(ff)
b_retimefig(time)

subplot('Position',[0.13 0.09 0.775 0.04])  %bar showing theta segments
fml = find(maxloc>pwind1&maxloc<pwind2);
if ~isempty(fml)
    draw(fml)
end
axis([0 x_lim(2) 0 1])
b_retimefig(time)

% Output matrix
OM = zeros(2,sw2);
OM(1,:) = maxloc;
for p = 1:sw2
    OM(2,p) = power(maxloc(p),p);
end

% ----------------------------------------------------------------------------------
% makes disjuct intervals from a set of points:
% points following each other continuusly form an interval, ie. all data points
% of the interval are in the set, and non of the data points out of the intervals
% are in the set
function draw(ip)
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
prepa = zeros(2,lenfdr+1);
prepa(1,1) = ip(1);
for t = 1:lenfdr
    prepa(2,t) = ip(fdr(t));
    prepa(1,t+1) = ip(fdr(t)+1);
end
prepa(2,end) = ip(end);
zs = zeros(1,lenfdr+1);
os = ones(1,lenfdr+1);
patcha = [prepa(1,:); prepa(1,:); prepa(2,:); prepa(2,:); prepa(1,:)];
patchb = [zs; os; os; zs; zs];
ptch = patch(patcha,patchb,'blue');
set(ptch,{'edgecolor'},get(ptch,{'facecolor'}));

% ----------------------------------------------------------------------------------
function b_retimefig(time)      %changes time labels to seconds
a0 = get(gca,'XTick');
iszero = 0;
if a0(1) == 0
    a0(1) = 1;
    iszero = 1;
end
a1 = time(a0);
if iszero
    a1(1) = 0;
end
set(gca,'XTickLabel',num2str(a1'))