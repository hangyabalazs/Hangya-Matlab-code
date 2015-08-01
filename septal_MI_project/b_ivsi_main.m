function [burstloc,intervals,minimums] = b_ivsi_main;
%IVSI_MAIN  Main function of IVSI.
%   IVSI_MAIN has three output arguments: first contains the localization of the first spikes of
%   bursts, second contains the localization of the first and last points of the theta intervals
%   of the eeg, third contains the localization of the eeg minimums.
%
%   See also IVSI.

% Input arguments check
error(nargchk(0,0,nargin));

% Loc. of first spikes of bursts
[burstloc] = b_burst_cls_for_ivsi3;

% Loc. of first and last points of theta intervals, loc. of minimums
[intervals,minimums] = b_wavephase_for_ivsi;

% Input variables
global IN
data = IN{1};
eeg = IN{2};
fname = IN{3};
pathname = IN{4};
datinx1 = IN{5};
datinx2 = IN{6};
time = IN{7};
unit = IN{8};
dt = IN{9};
meret = IN{10};
mintafr = IN{11};
xlimit = IN{12};

% Ivsi plot
ooo = figure;
lengthint = zeros(1,length(burstloc));
lengthint_b = zeros(1,length(burstloc));
lengthburst = zeros(1,length(burstloc));
lengthburst_b = zeros(1,length(burstloc));
for i = 2:length(burstloc), %burst interval begins before the theta interval (red)
    c = burstloc(i);
    c1 = find(intervals<c);
    if rem(length(c1),2)==1,
        m1 = find(minimums<c);
        if isempty(m1) == 0,
            m2 = m1(end);
            if m2 < length(minimums),
                d1 = find(intervals<minimums(m2));
                d2 = find(intervals<minimums(m2+1));
                if length(d1) == length(d2),
                    lengthint_b(i) = (minimums(m2+1) - minimums(m2))/10;
                    lengthburst_b(i) = (burstloc(i) - burstloc(i-1))/10;
                    plot3(lengthint_b(i),lengthburst_b(i),burstloc(i),'rs');
                    hold on;
                end;
            end;
        end;
    end;
end;
o = figure;
plot(lengthint_b,lengthburst_b,'rs');
hold on;
set(0,'CurrentFigure',ooo);
for i = 1:length(burstloc)-1, %burst interval begins after the theta interval (blue)
    c = burstloc(i);
    c1 = find(intervals<c);
    if rem(length(c1),2)==1,
        m1 = find(minimums<c);
        if isempty(m1) == 0,
            m2 = m1(end);
            if m2 < length(minimums),
                d1 = find(intervals<minimums(m2));
                d2 = find(intervals<minimums(m2+1));
                if length(d1) == length(d2),
                    lengthint(i) = (minimums(m2+1) - minimums(m2))/10;
                    lengthburst(i) = (burstloc(i+1) - burstloc(i))/10;
                    plot3(lengthint(i),lengthburst(i),burstloc(i),'bo');
                    for i = 1:length(lengthburst)-1,   %horizontal lines
                        if lengthburst(i) ~= 0 & lengthburst_b(i+1) ~=0,
                            line([lengthint(i) lengthint_b(i+1)],...
                                [lengthburst(i) lengthburst_b(i+1)],...
                                [burstloc(i) burstloc(i+1)]);
                        end;
                    end;
                    for i = 1:length(lengthburst),   %vertical lines
                        if lengthburst(i) ~= 0 & lengthburst_b(i) ~=0,
                            line([lengthint(i) lengthint_b(i)],...
                                [lengthburst(i) lengthburst_b(i)],...
                                [burstloc(i) burstloc(i)]);
                        end;
                    end;
                    hold on;
                end;
            end;
        end;
    end;
end;
grid on;
set(0,'CurrentFigure',o);
plot(lengthint,lengthburst,'bo');
hold on;
for i = 2:length(lengthburst),
    if lengthburst(i-1) == 0,
        plot(lengthint(i),lengthburst(i),'b+');
    end;
end;
for i = 2:length(lengthburst_b),
    if lengthburst_b(i-1) == 0,
        plot(lengthint_b(i),lengthburst_b(i),'rx');
    end;
end;
axis([0 800 0 800]);
line([0 800],[0 800],'Color','g');
fname2 = [fname(1:3) ' ' fname(5:6)];
title(['1st spikes vs theta intervals ' fname2  ' ' num2str(datinx1/10000) ...
        ':' num2str(datinx2/10000)],'Color','red');
xlabel('theta intervals');
ylabel('time intervals between 1st spikes of bursts');

% Ivsi plot with lines connecting every interval with the following
oo = copyobj(o,0);
for i = 1:length(lengthburst)-1,   %horizontal lines
    if lengthburst(i) ~= 0 & lengthburst_b(i+1) ~=0,
        line([lengthint(i) lengthint_b(i+1)], [lengthburst(i) lengthburst_b(i+1)]);
    end;
end;
for i = 1:length(lengthburst),   %vertical lines
    if lengthburst(i) ~= 0 & lengthburst_b(i) ~=0,
        line([lengthint(i) lengthint_b(i)], [lengthburst(i) lengthburst_b(i)]);
    end;
end;

% Plotting the two interval variable separated
figure;
lengthint2 = lengthint(find(lengthint));
lengthburst2 = lengthburst(find(lengthburst));
ti = [1:length(lengthint2)];
plot(ti,lengthint2,'-ro',ti,lengthburst2,'-bx');    %red is theta

% Square Signal Plot of the two interval variables: the two variables are ploted separated,
% where hight corresponds with the value
figure;
du1 = zeros(1,length(time));
for i = 1:length(minimums)-1,
    du1(minimums(i):minimums(i+1)) = minimums(i+1) - minimums(i); 
end;
du2 = zeros(1,length(time));
for i=1:length(burstloc)-1,
    du2(burstloc(i):burstloc(i+1)) = burstloc(i+1) - burstloc(i); 
end;
plot(du1,'-r'); %red is theta
hold on;
plot(du2,'-b');