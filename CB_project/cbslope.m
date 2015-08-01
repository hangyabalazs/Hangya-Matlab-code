function [slope amp] = cbslope(data)
%CBSLOPE   Calculates slopes and amplitudes for stimulus response.
%   [S, A] = CBSLOPE(DATA) calculates slopes (S) and amplitudes (A) for 
%   stimulus response in DATA. DATA should contain a channel of stimuli.
%   Amplitude is defined as the difference of prestimulus data value and
%   response minimum. A line is fitted on the descending phase of the
%   response where the L2-error of the fitting is minimal when lines of one
%   third length of the descnding phase are used. Slope of this line is 
%   returned in the output argument.   
%
%   See also LINEFIT.

% Extract data
eeg = data(:,1)';
stim = data(:,2)';
sr = 20000;

% Filtering
% nqf = sr / 2;
% flt = fir1(512,10/nqf,'high');
% feeg = filtfilt(flt,1,eeg);
% 
% figure;plot(feeg)
% axis([180000 200000 -1.1 0.5])

% Fit line
vdisc = find(diff(stim)>0.9);
lv = length(vdisc);
amp = zeros(1,lv);
slope = zeros(1,lv);
intercept = zeros(1,lv);
for k = 1:lv
    sp = vdisc(k);
    wn = eeg(sp+0.001*sr:round(sp+0.025*sr)); %modified from 0.0005 to 0.001
    mn = min(wn);
    mninx = find(wn==mn);
    mninx = mninx(end);
    mx = max(wn(1:mninx));
    mxinx = find(wn(1:mninx)==mx);
    mxinx = mxinx(1);
    linlen = round((mninx-mxinx)/3);
    time = linspace(0,linlen/sr,linlen+1);
    lm = mninx-mxinx-linlen;
    r = zeros(1,lm);
    gr = zeros(1,lm);
    icp = zeros(1,lm);
    for m = 1:lm
        [gr0,icp0,err0] = linefit(time,wn(mxinx+m:mxinx+m+linlen));
        r(m) = err0;
        gr(m) = gr0;
        icp(m) = icp0;
    end
    lsr = min(r);
    grinx = find(r==lsr);
    slope(k) = gr(grinx) / 1000;   % slope is given in mV/ms
    intercept(k) = icp(grinx);
    amp(k) = eeg(sp) - mn;
    
%     figure
%     plot(wn)
%     hold on
%     t = [mxinx+grinx:mxinx+grinx+linlen];
%     y = time .* slope(k) + intercept(k);
%     plot(t,y,'r')
end