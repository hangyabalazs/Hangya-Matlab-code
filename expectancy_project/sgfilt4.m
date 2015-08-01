function sgfilt4(data)
%SGFILT4    Filter test.
%   Filter test for expectancy data.
%
%   See also SGFILT, SGFILT2 and SGFILT3.

% Load data
global DATAPATH
resdir = [DATAPATH 'Expectancy\filter2\Fir1\'];
cd(resdir)

% Filtering
dim2 = 4;       % 91%
dim3 = 2;       % Cz
for dim1 = 1:13
    dat = squeeze(data(dim1,dim2,dim3,:,:))';
    sr = 1000;
    nqf = sr / 2;
    fraw = zeros(100,3000);
    nraw = zeros(100,3000);
    ahraw = zeros(100,3000);
    flt = fir1(512,[0.5 3]/nqf);
    for k = 1:100
        raw = dat(k,:);
        raw(1:1550) = randn(1,1550) + raw(1550);        % first 1550 points replaced with random noise
        nraw(k,:) = raw;
        fraw(k,:) = filtfilt(flt,1,raw);
        ahraw(k,:) = angle(hilbert(fraw(k,:)));
    end
    frr = mean(fraw);
    H = figure;
    plot(mean(nraw))
    hold on
    plot(frr,'r')
    saveas(H,['ksz' num2str(dim1) '.fig'])
end