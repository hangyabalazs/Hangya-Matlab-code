function b_ward_clustering
%WARD_CLUSTERING    Ward,s cluster analysis for unit.
%   WARD_CLUSTERING is a version of BURST_CLS. It saves a histogram of the intraburst instant
%   frequency and saves also two matrices: burst length in 'bl' and intraburst spike number in
%   'sn'. Output is stored in d:\matlab6_1\data\burst_clusters\ufos2 directory.
%
%   The program calls a subfunction called PLOTBURST who saves the 'burst plot' in the current
%   directory. Note that you should change the current directory if you don't want to see the
%   plots in the directory you are in while running WARD_CLUSTERING itself (not by CLUSTERRUN).
%
%   See also BURST_CLS and CLUSTERRUN.

% Input arguments check
error(nargchk(0,0,nargin));

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

% Discrimination variables
global DISC
output = DISC{2};
vdisc = DISC{3};
kuszob = DISC{4};
instfrek = DISC{5};
isi = DISC{6};

% Return plot
ivs=diff(vdisc);
livs=length(ivs);
d1=ivs(1:livs-1);
d2=ivs(2:livs);
figure
plot(d1,d2,'.')
title('return plot of ISIs')
xlabel('nth ISI')
ylabel('(n+1)th ISI')
k=waitforbuttonpress;
if k==0 close; end;

% Ward's clustering
dmtx=zeros(livs,livs);
iivs=ivs';
dmtx(1,:)=ivs;
dmtx(:,1)=iivs;
dist=pdist2(dmtx);
links=linkage(dist,'ward');
[q w]=dendrogram(links);
c3=cluster(links,3);
c4=cluster(links,4);
c5=cluster(links,5);
c6=cluster(links,6);
c7=cluster(links,7);
c31=find(c3==1);
c32=find(c3==2);
c41=find(c4==3);
c42=find(c4==2);
c43=find(c4==1);
c51=find(c5==4);
c52=find(c5==3);
c53=find(c5==2);
c54=find(c5==1);
c61=find(c6==5);
c62=find(c6==4);
c63=find(c6==3);
c64=find(c6==2);
c65=find(c6==1);
c71=find(c7==6);
c72=find(c7==5);
c73=find(c7==4);
c74=find(c7==3);
c75=find(c7==2);
c76=find(c7==1);
a1=length(c31); a2=length(c32);
a3=length(c41); a4=length(c42); a5=length(c43);
a6=length(c51); a7=length(c52); a8=length(c53); a9=length(c54);
a10=length(c61); a11=length(c62); a12=length(c63); a13=length(c64); a14=length(c65);
a15=length(c71); a16=length(c72); a17=length(c73); a18=length(c74); a19=length(c75); a20=length(c76);
alla3=[a1 a2];
alla4=[a3 a4 a5];
alla5=[a6 a7 a8 a9];
alla6=[a10 a11 a12 a13 a14];
alla7=[a15 a16 a17 a18 a19 a20];
maxa3=find(alla3==max(alla3));
maxa4=find(alla4==max(alla4));
maxa5=find(alla5==max(alla5));
maxa6=find(alla6==max(alla6));
maxa7=find(alla7==max(alla7));
binum=fix(0.606+0.4*(livs-1));
if maxa3==1 ivss3=ivs(c31); restivs3=ivs(c32);
else ivss3=ivs(c32); restivs3=ivs(c31); end;
binum1=fix(0.606+0.4*(length(ivss3)-1));
[n1 m1]=hist(ivss3,binum1);
if maxa4==1 ivss4=ivs(c41); restivs4=ivs([c42' c43']); end;
if maxa4==2 ivss4=ivs(c42); restivs4=ivs([c41' c43']); end;
if maxa4==3 ivss4=ivs(c43); restivs4=ivs([c41' c42']); end;
binum2=fix(0.606+0.4*(length(ivss4)-1));
[n2 m2]=hist(ivss4,binum2);
if maxa5==1 ivss5=ivs(c51); restivs5=ivs([c52' c53' c54']); end;
if maxa5==2 ivss5=ivs(c52); restivs5=ivs([c51' c53' c54']); end;
if maxa5==3 ivss5=ivs(c53); restivs5=ivs([c51' c52' c54']); end;
if maxa5==4 ivss5=ivs(c54); restivs5=ivs([c51' c52' c53']); end;
binum3=fix(0.606+0.4*(length(ivss5)-1));
[n3 m3]=hist(ivss5,binum3);
if maxa6==1 ivss6=ivs(c61); restivs6=ivs([c62' c63' c64' c65']); end;
if maxa6==2 ivss6=ivs(c62); restivs6=ivs([c63' c64' c65' c61']); end;
if maxa6==3 ivss6=ivs(c63); restivs6=ivs([c61' c62' c64' c65']); end;
if maxa6==4 ivss6=ivs(c64); restivs6=ivs([c61' c62' c63' c65']); end;
if maxa6==5 ivss6=ivs(c65); restivs6=ivs([c61' c62' c63' c64']); end;
binum4=fix(0.606+0.4*(length(ivss6)-1));
[n4 m4]=hist(ivss6,binum4);
if maxa7==1 ivss7=ivs(c71); restivs7=ivs([c72' c73' c74' c75' c76']); end;
if maxa7==2 ivss7=ivs(c72); restivs7=ivs([c71' c73' c74' c75' c76']); end;
if maxa7==3 ivss7=ivs(c73); restivs7=ivs([c71' c72' c74' c75' c76']); end;
if maxa7==4 ivss7=ivs(c74); restivs7=ivs([c71' c72' c73' c75' c76']); end;
if maxa7==5 ivss7=ivs(c75); restivs7=ivs([c71' c72' c73' c74' c76']); end,
if maxa7==6 ivss7=ivs(c76); restivs7=ivs([c71' c72' c73' c74' c75']); end;
binum5=fix(0.606+0.4*(length(ivss7)-1));
[n5 m5]=hist(ivss7,binum5);
[ng mg]=hist(ivs,binum);
lmg=fix(max(mg));
figure
subplot(2,1,1)
bar(mg,ng)
subplot(2,1,2)
bar(m1,n1)
set(gca,'XLim',[0 lmg])
figure
subplot(2,1,1)
bar(mg,ng)
subplot(2,1,2)
bar(m2,n2)
set(gca,'XLim',[0 lmg])
figure
subplot(2,1,1)
bar(mg,ng)
subplot(2,1,2)
bar(m3,n3)
set(gca,'XLim',[0 lmg])
figure
subplot(2,1,1)
bar(mg,ng)
subplot(2,1,2)
bar(m4,n4)
set(gca,'XLim',[0 lmg])
figure
subplot(2,1,1)
bar(mg,ng)
subplot(2,1,2)
bar(m5,n5)
set(gca,'XLim',[0 lmg])
dec=input('Choose the No of clusters you like best -or 0 for none:-)');
if dec==0
%     clear c31 c41 c51 c61 c71 ivss3 ivss4 ivss5 ivss6 ivss7
    disp('Your cell is not bursty -such a shame:-(')
    bursting=0;
    burst=[];
    plotburst(bursting,burst)
    return    
end;
bins=[10:10:500];
where='d:\matlab6_1\data\burst_clusters\ufos2\';
segm_ident=['_' num2str(datinx1/10000) '_' num2str(datinx2/10000)];
fnam=[where fname(1:6) segm_ident];
if dec==3
%     clear c41 c51 c61 c71 ivss4 ivss5 ivss6 ivss7;
    bursting=1;
    blf=1/(max(ivss3)/10000);
    aveintra=1/(mean(ivss3)/10000);
    stdintra=1/(std(ivss3)/100000);
    cov_burst=std(ivss3)/mean(ivss3);
    cov_non_burst=std(restivs3)/mean(restivs3);
    ibfs3=(1./ivss3)*10000;
    ibfs3=ibfs3';
    [f3 g3]=hist(ibfs3,bins);
    fsave3=figure; bar(g3,f3); saveas(fsave3,fnam,'fig');
    save(fnam,'ibfs3','-ASCII');
end;
if dec==4
%     clear c31 c51 c61 c71 ivss3 ivss5 ivss6 ivss7;
    bursting=1;
    blf=1/(max(ivss4)/10000);
    aveintra=1/(mean(ivss4)/10000);
    stdintra=1/(std(ivss4)/100000);
    cov_burst=std(ivss4)/mean(ivss4);
    cov_non_burst=std(restivs4)/mean(restivs4);
    ibfs4=(1./ivss4)*10000;
    ibfs4=ibfs4';
    [f4 g4]=hist(ibfs4,bins);
    fsave4=figure; bar(g4,f4); saveas(fsave4,fnam,'fig');
    save(fnam,'ibfs4','-ASCII');
end;
if dec==5
%     clear c31 c41 c61 c71 ivss3 ivss4 ivss6 ivss7;
    bursting=1;
    blf=1/(max(ivss5)/10000);
    aveintra=1/(mean(ivss5)/10000);
    stdintra=1/(std(ivss5)/100000);
    cov_burst=std(ivss5)/mean(ivss5);
    cov_non_burst=std(restivs5)/mean(restivs5);
    ibfs5=(1./ivss5)*10000;
    ibfs5=ibfs5';
    [f5 g5]=hist(ibfs5,bins);
    fsave5=figure; bar(g5,f5); saveas(fsave5,fnam,'fig');
    save(fnam,'ibfs5','-ASCII');
end;
if dec==6
%     clear c31 c41 c51 71 ivss3 ivss4 ivss5 ivss7;
    bursting=1;
    blf=1/(max(ivss6)/10000);
    aveintra=1/(mean(ivss6)/10000);
    stdintra=1/(std(ivss6)/100000);
    cov_burst=std(ivss6)/mean(ivss6);
    cov_non_burst=std(restivs6)/mean(restivs6);
    ibfs6=(1./ivss3)*10000;
    ibfs6=ibfs6';
    [f6 g6]=hist(ibfs6,bins);
    fsave6=figure; bar(g6,f6); saveas(fsave6,fnam,'fig');
    save(fnam,'ibfs6','-ASCII');
end;
if dec==7
%     clear c31 c41 c51 c61 ivss3 ivss4 ivss5 ivss6;
    bursting=1;
    blf=1/(max(ivss7)/10000);
    aveintra=1/(mean(ivss7)/10000);
    stdintra=1/(std(ivss7)/100000);
    cov_burst=std(ivss7)/mean(ivss7);
    cov_non_burst=std(restivs7)/mean(restivs7);
    ibfs7=(1./ivss7)*10000;
    ibfs7=ibfs7';
    [f7 g7]=hist(ibfs7,bins);
    fsave7=figure; bar(g7,f7); saveas(fsave7,fnam,'fig');
    save(fnam,'ibfs7','-ASCII');
end;
clear alla* maxa*

burstout=struct('blf',{0},'ibsn',{0},'ibfr',{0},'blen',{0}, ...
                'aibsn',{0},'ablen',{0}, ...
                'sibfr',{0},'sibsn',{0},'sblen',{0});

% Marking of the bursts & burst characteristics
if bursting,
liv=find(ivs<10000/blf);
liv1=[-1 liv]; liv=[liv 0];
bliv=find(liv~=liv1+1);
blivv=bliv(2:end)-1; bliv(end)=[];
burst=[liv(bliv);liv(blivv)+1];
fsp=vdisc(burst(1,:));               %1st spikes of bursts
bn=size(burst,2);
burstout.blf=blf;
ibsn=burst(2,:)+1-burst(1,:); ibsn=ibsn';
bl=(vdisc(burst(2,:))-vdisc(burst(1,:)))/10; bl=bl';
fnam2=[fnam '_bl'];
fnam3=[fnam '_sn'];
save(fnam2,'bl','-ASCII');
save(fnam3,'ibsn','-ASCII');

burstout.ibfr=aveintra;
burstiness=sum(burstout.ibsn)/length(vdisc);
burstout.ablen=sum(burstout.blen)/(10*bn);
burstout.aibsn=sum(burstout.ibsn)/bn;
burstout.sblen=std(burstout.blen)/(10*bn);
burstout.sibfr=stdintra;
burstout.sibsn=std(burstout.ibsn)/bn;
bfs=burst(1,:);
blength=burstout.blen;
[t1 t2]=hist(blength,15);
end;
ibn=mean(burstout.ibsn);
disp(ibn);
disp(burstout.ablen);
disp(burstiness);
disp(burstout.ibfr);
disp(cov_burst);
disp(cov_non_burst);
dt=0.0001;
time=[0:length(unit)-1]*dt;
close all
plotburst(bursting,burst)

% ---------------------------------------------------------------------------------
function plotburst(bursting,burst)

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

% Discrimination variables
global DISC
output = DISC{2};
vdisc = DISC{3};
kuszob = DISC{4};
instfrek = DISC{5};
isi = DISC{6};

% Plot and save
hangya = figure;
plot(time,unit,'m');
if bursting==1
    hold on;
    mi=kuszob; ma=kuszob+0.5;
    for j=1:length(burst)
        rajz=[time(vdisc(burst(1,j))) time(vdisc(burst(2,j)));mi ma];
        line(rajz(1,:),rajz(2,:),'color','b');
    end;
    hold off;
end;
title('Bursts');
filenam = fname(1:6);
eval(['saveas(hangya,''BURSTS_',filenam,'.fig'')']);