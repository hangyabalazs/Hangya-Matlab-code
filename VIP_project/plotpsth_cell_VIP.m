%%

PSTHs = allpsth(NTactinx,:);
NumPsths = size(PSTHs,1);

%%


for k =  1:NumPsths
    PSTHs(k,:) = PSTHs(k,:) - mean(PSTHs(k,:));
end

%%

mxinx = nan(1,NumPsths);
for k = 1:NumPsths
    mxinx(k) = sum(PSTHs(k,:)>0);
end

%%

mxinx = nan(1,NumPsths);
fr1 = nan(1,NumPsths);
fr2 = nan(1,NumPsths);
tinx1 = time > -400 & time < -150;
tinx2 = time > 0 & time < 300;
for k = 1:NumPsths
    fr1(k) = sum(PSTHs(k,tinx1));
    fr2(k) = sum(PSTHs(k,tinx2));
    mxinx(k) = fr2(k) - fr1(k);
end

%%

mxinx = nan(1,NumPsths);
tinx = time > 20 & time < 30;
for k = 1:NumPsths
    mxinx(k) = sum(PSTHs(k,tinx));
end

%%

[mx mxinx] = max(PSTHs,[],2);
[srt srtinx] = sort(mx,'descend');
figure
imagesc(time,1:NumPsths,PSTHs(srtinx,:))

%%

figure
pinx = time > -300 & time < 300;
imagesc(time(pinx),1:NumPsths,PSTHs(srtinx,pinx))

%%

fr1 = nan(1,NumPsths);
fr2 = nan(1,NumPsths);
tinx1 = time > 0 & time < 100;
tinx2 = time > 100 & time < 200;
for k = 1:NumPsths
    fr1(k) = sum(PSTHs(k,tinx1));
    fr2(k) = sum(PSTHs(k,tinx2));
end

%%

fr1 = nan(1,NumPsths);
fr2 = nan(1,NumPsths);
tinx1 = time > 0 & time < 100;
tinx2 = time > 100 & time < 200;
for k = 1:NumPsths
    [fr1(k) mxl] = max(abs(PSTHs(k,tinx1)));
    if PSTHs(k,mxl) < 0
        fr1(k) = -fr1(k);
    end
    [fr2(k) mxl] = max(abs(PSTHs(k,tinx2)));
    if PSTHs(k,mxl) < 0
        fr2(k) = -fr2(k);
    end
end

%%

figure;
plot(fr1,fr2,'.')

%%

figure;
plot(fr1+fr2,fr1-fr2,'.')

%%

fr1 = nan(1,NumPsths);
fr2 = nan(1,NumPsths);
tinx1 = time >= 20 & time <= 32;
tinx2 = time > 32 & time < 150;
for k = 1:NumPsths
    fr1(k) = mean(PSTHs(k,tinx1));
    fr2(k) = mean(PSTHs(k,tinx2));
end