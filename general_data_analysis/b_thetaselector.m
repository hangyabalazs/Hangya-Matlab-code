function H = b_thetaselector(wave,f)
%THETASELECTOR   Theta selection upon wavelet power maximum localization.
%   H = THETASELECTOR2(WAVE,F) find maximum localization of WAVE power in each 
%   coloumn, using F scalevector. Only maximums falling into the highest 15 % 
%   of the power values are considered ('real maxes').
%   WAVE power is plotted on figure H with white pots signaling the real maximums.
%   A bar shows the segments with real maximums, and another one shows the theta
%   segments (segments where real maximums fall in the theta band).
%
%   Distribution of theta power and all power values are plotted in another figure
%   window.
%
%   Theta band: 2.5 - 6 Hz.
%
%   See also THETASELECTOR3.

% Input argument check
error(nargchk(2,2,nargin));

% Finding theta band
fnd = find(f>6);
pwind1 = fnd(end);
fnd = find(f<3);
pwind2 = fnd(1);

% Distribution of theta power
power = (abs(wave)) .^ 2;
thetapower = power(pwind1:pwind2,:);

lentp = size(thetapower,1) * size(thetapower,2);
bno = 100;
figure
thetapower2 = reshape(thetapower,1,lentp);
[x y] = hist(thetapower2,bno);
bar(y,x);
title('Distribution of theta power');

% Distribution of all power
sw1 = size(power,1);
sw2 = size(power,2);
lenw = sw1 * sw2;
bno = 100;
figure
power2 = reshape(power,1,lenw);
[x y] = hist(power2,bno);
bar(y,x);
title('Distribution of all power');

% Maximum localizations & real maximums
maxes = max(power);
maxloc = zeros(1,sw2);
for i = 1:sw2
    maxloc(i) = find(power(:,i)==maxes(i));
end
bno = [1:sw1];
H = figure;
[x y] = hist(maxloc,bno);
bar(y,x);

spower = sort(power2);
percent = fix(0.15 * lenw);  % selecting the highest 15 per cent
realmaxes = spower(end-percent+1:end);
realmaxloc = find(maxes>=realmaxes(1));

% Plot
figure
subplot('Position',[0.13 0.11+0.15 0.775 0.815-0.15])
imagesc(power)
hold on
maxloc2 = zeros(size(maxloc));
maxloc2(realmaxloc) = maxloc(realmaxloc);
plot(maxloc2,'Color',[1 1 1],'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 1 1],...
    'MarkerSize',0.5,'LineStyle','none')
x_lim = get(gca,'XLim');
line([x_lim(1) x_lim(2)],[pwind1 pwind1],'Color','g');
line([x_lim(1) x_lim(2)],[pwind2 pwind2],'Color','g');

subplot('Position',[0.13 0.17 0.775 0.04])  %bar showing segments with real maxes
draw(realmaxloc)
axis([0 x_lim(2) 0 1])

subplot('Position',[0.13 0.08 0.775 0.04])  %bar showing theta segments
frml = find(maxloc(realmaxloc)>pwind1&maxloc(realmaxloc)<pwind2);
rfmrl = realmaxloc(frml);
if ~isempty(rfmrl)
    draw(rfmrl)
end
axis([0 x_lim(2) 0 1])

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