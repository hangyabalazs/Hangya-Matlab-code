%% plot MI

aMI = mean(MI,3);
figure;
imagesc(aMI)
b_rescaleaxis('X',f(1:5:end))
b_rescaleaxis('Y',f(1:5:end))