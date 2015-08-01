function [H,slope] = somslope(DataPoints)

% Slice
n = length(DataPoints);



wn = DataPoints;
mxinx = find(wn==max(wn));
mxinx = mxinx(end);
linlen = 50;
lm = min(n-linlen,mxinx-linlen+1);
lsrate = 1000;
time = linspace(0,linlen/lsrate,linlen+1);
r = zeros(1,lm);
gr = zeros(1,lm);
icp = zeros(1,lm);
for m = 1:lm
    [gr0,icp0,err0] = ...
        linefit(time*1000,wn(m:m+linlen));    % time is rescaled because large scaling difference screws up linefit
    r(m) = err0;
    gr(m) = gr0 * 1000;    % rescale to agree with original time vector
    icp(m) = icp0;
end
r(gr<0) = -1000;   % looking for rising edge, therefore set error of descending segments to -inf
grinx = find(gr==max(gr));
if isempty(grinx);
    slope = 500;
    intercept = 0;
else
    grinx = grinx(end);
    slope = gr(grinx);
    intercept = icp(grinx);
end

H = [];
H = figure;
plot(wn)
hold on
t = grinx:grinx+linlen;
y = time .* slope + intercept;
if ~isempty(grinx)
    plot(t,y,'r')
    r(grinx)
end