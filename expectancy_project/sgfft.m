%% STFT

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
y = zeros(100,501);
w = zeros(100,501);
for k = 1:100
    raw = dat(k,1000:1500);
    [y(k,:),w(k,:)] = b_fft2(raw,sr);
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
W = mean(w);
figure
plot(W(2:25),Y(2:25))

%% STFT

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
us = 4;
nqf = sr / 2;
y = zeros(100,501*us);
w = zeros(100,501*us);
for k = 1:100
    raw = dat(k,1000:1500);
    raw2 = resample(double(raw),us,1);
    [y(k,:),w(k,:)] = b_fft2(raw2,sr*us);
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
W = mean(w);
figure
plot(W(2:100),Y(2:100))

%% STFT

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
rp = 4;
y = zeros(100,501*rp);
w = zeros(100,501*rp);
for k = 1:100
    raw = dat(k,1000:1500);
    raw2 = repmat(raw,1,rp);
    [y(k,:),w(k,:)] = b_fft2(raw2,sr);
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
W = mean(w);
figure
plot(W(2:100),Y(2:100))

%% STFT - ultimate

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
nfft = 2004;
klim1 = 1;
klim2 = 100;
knum = klim2-klim1+1;
y = zeros(knum,nfft);
for k = klim1:klim2
    raw = dat(k,1000:1500);
    [y(k,:),w] = b_fft2(raw,sr,nfft);
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
figure
plot(w(2:100),Y(2:100))
% figure
% imagesc(y(:,5:40),[0 10^6])
% figure
% mesh(y(:,5:40))

%% STFT

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
nfft = 2004;
y = zeros(100,nfft);
w = zeros(100,nfft);
for k = 1:100
    raw = dat(k,1000:1500);
    [y(k,:)] = abs(fft(raw2,nfft));
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
figure
plot(Y(2:1000))

%% STFT

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
y = zeros(100,1001);
w = zeros(100,1001);
for k = 1:100
    raw = dat(k,500:1500);
    [y(k,:),w(k,:)] = b_fft2(raw,sr);
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
W = mean(w);
figure
plot(W(2:100),Y(2:100))


%% STFT - 1-500 ms

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
nfft = 2004;
klim1 = 1;
klim2 = 100;
knum = klim2-klim1+1;
y = zeros(knum,nfft);
for k = klim1:klim2
    raw = dat(k,1:501);
    [y(k,:),w] = b_fft2(raw,sr,nfft);
%     figure
%     plot(w(2:25),y(2:25))
end
Y = mean(y);
figure
plot(w(2:100),Y(2:100))
% figure
% imagesc(y(:,5:40),[0 10^6])
% figure
% mesh(y(:,5:40))

%% STFT - mean data fft

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(erpdata(13,dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
nfft = 4000;
klim1 = 1;
klim2 = 100;
knum = klim2-klim1+1;
raw = mean(dat(klim1:klim2,1000:1500));
[y,w] = b_fft2(raw,sr,nfft);
figure
plot(w(2:100),y(2:100))
figure
plot(w(2:100),log(y(2:100)))

%% STFT - ultimate

dim1 = 4;       % condition
dim2 = 3;       % Cz
dat = squeeze(erpdata(:,dim1,dim2,:,:));
sr = 1000;
nqf = sr / 2;
nfft = 4004;
klim1 = 1;
klim2 = 100;
snum = 13;
knum = klim2-klim1+1;
y = zeros(snum,knum,nfft);
for subj = 1:snum;
    for k = klim1:klim2
        raw = dat(subj,1000:1500,k);
        [y(subj,k,:),w] = b_fft2(raw,sr,nfft);
        %     figure
        %     plot(w(2:25),y(2:25))
    end
    figure;
    plot(w(1:100),squeeze(mean(y(subj,:,1:100),2)))
    title(num2str(subj))
end
% Y = squeeze(mean(y,2));
% Y = squeeze(Y(6,:));
figure
y2 = squeeze(reshape(y,1,1300,4004));
imagesc(w(1:100),1:1300,y2(:,1:100),[0 10^7])
Y = squeeze(mean(mean(y,2),1));
figure
plot(w(1:100),log(Y(1:100)))
figure
plot(w(2:100),Y(2:100))

%% STFT - subjects

dim2 = 3;       % Cz
sr = 1000;
nqf = sr / 2;
nfft = 4004;
klim1 = 1;
klim2 = 100;
snum = 13;
knum = klim2-klim1+1;
y = zeros(snum,knum,nfft);
clr = [0 0 1;[255 153 51]/256;[0 1 0];[1 0 1]];
for dim1 = 1:4
    dat = squeeze(erpdata(:,dim1,dim2,:,:));
    for subj = 1:snum;
        for k = klim1:klim2
            raw = dat(subj,1000:1500,k);
            [y(subj,k,:),w] = b_fft2(raw,sr,nfft);
            %     figure
            %     plot(w(2:25),y(2:25))
        end
        figure(subj);
        hold on
        plot(w(1:100),squeeze(mean(y(subj,:,1:100),2)),'Color',clr(dim1,:),'LineWidth',2)
        title(num2str(subj))
    end
end