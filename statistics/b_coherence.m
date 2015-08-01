%COHERENCE   Imported from Viktor.

% sdf_cohere
% if ~exist('meghivott'),
%     in;
%     disc_futtat;
% end;

meghivott=1;
sdf_f;
coh=cohere(spdenz,eeg,16384,10000);
tagfrar=16384/10000;
tag20=fix(20*tagfrar);
xax=[0:tag20]/tagfrar;
tag3=min(find(xax>2.5)); tag5=max(find(xax<5.5));
fis.coh=figure;
    plot(xax,coh(1:tag20+1));
set(gca,'ylim',[0 1]);
xlabel('frequency (Hz)');
hold on;
for i=tag3:tag5,
    rajz=[xax(i) xax(i);0 coh(i)];
    line(rajz(1,:),rajz(2,:),'color','r');
end;
hold off;
        
[maxcoh mcfr]=max(coh(tag3:tag5));
mcfr=xax(mcfr+tag3-1);

clear meghivott;

% finame=[figmentpath fname(1:6) '_coh_' iv];
% saveas(fis.coh,finame,'tiffn');
% fprintf(fresult,[iv '_COH_max_freq_in_theta_range %11e\n'],mcfr);
% fprintf(fresult,[iv '_COH_max_freq_value          %11e\n'],maxcoh);
        