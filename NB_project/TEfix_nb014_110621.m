son3=son2-son2(1)+ts(1);
figure;plot(ts);hold on;plot(son3,'r')
ts2=ts;
ts(30:38)=[];
figure;plot(ts);hold on;plot(son3,'r')

tso=TE2.StimulusOn;
tso(30:38)=[];
TE2.TrialStart = son2 - tso;

fnm=fieldnames(TE2);
for k=1:length(fieldnames(TE2))
    if length(TE2.(fnm{k}))==289
        TE2.(fnm{k})(30:38)=[];
    end
end

save([sessionpath filesep 'TrialEvents.mat'],'-struct','TE2')