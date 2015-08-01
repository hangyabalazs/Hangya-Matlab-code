function b_burst_cls_auto
%BURST_CLS_AUTO  Automatic version of BURST_CLS.
%   The program uses five criteria to decide whether the cell shows thetha frequency
%   bursting. First of all the spike train should contain more than ten spikes. Second, 
%   relative theta power should be over 25 per cent. Third, the shortest interval 
%   should be in the largest cluster. Fourth, no interval longer than 150 ms should 
%   be in that cluster. Fifth, none of the bursts should be longer than 300 ms.
%
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code.
%
%   See also BURST_CLS.

% Input arguments check
error(nargchk(0,0,nargin)); 

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'analysenow5\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
cd([DATAPATH,'relfreq_figs_proba\']);  %Here are the results
npt = input('Discard existing content or append data while writing burstauto.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('burstauto.txt','w');
elseif npt == 'a',
    fid = fopen('burstauto.txt','a');
else
    error('Unexpected answer for input.')
end;

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    meret=size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2)==1,
        data=[data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
   
% Searching for the matching datinx file
    cl = fname(1:6);
    filename = [cl,'.txt'];
    pont = findstr(filename,'.');
    filenam = filename(1:pont(1)-1);
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
    ww = 0;
    
% Computing the input variables
    for p = 1:3,
        datinx1 = v(2*p-1); %first point of the interval
        datinx2 = v(2*p); %last point of the interval
   
        dt = 0.0001;
        mintafr = 1 / dt;
        if datinx1(1) ~= 0 & datinx2(1) ~= 0,
            ww = ww + 1;
            unit = data(datinx1:datinx2,2);
            unit = unit';
            eeg = data(datinx1:datinx2,1);
            eeg = eeg';
            time = [0:length(unit)-1] * dt; 
            xlimit = [min(time),max(time)];
            global IN
            IN = cell(1,12);
            IN{1} = data;
            IN{2} = eeg;
            IN{3} = fname;
            IN{4} = pathname;
            IN{5} = datinx1;
            IN{6} = datinx2;
            IN{7} = time;
            IN{8} = unit;
            IN{9} = dt;
            IN{10} = meret;
            IN{11} = mintafr;
            IN{12} = xlimit;
            
% Discrimination                
            kuszob = v(6+ww);
            b_disc(kuszob);
            global DISC
            id = DISC{1};
            output = DISC{2};
            vdisc = DISC{3};
            kuszob = DISC{4};
            instfrek = DISC{5};
            isi = DISC{6};
            
% Spike Density Function
            ipunit=zeros(1,length(unit));
            ipunit(vdisc)=1;
            wbh=blackmanharris(1000);
            wipunit=conv(ipunit,wbh);
            wipunit=wipunit(1:length(unit));
            [punit funit]=periodogram(wipunit,blackmanharris(length(wipunit)),65536,10000);
            relthetapower = (sum(punit(21:41))) / (sum(punit(4:end)));    %relative theta power
            maxamp = max(punit(21:41));    %maximum theta power
            maxloc = find(punit==maxamp);
            maxfreq = funit(maxloc);    %the frequency where theta power is maximal
            pmin = min(wipunit);
            pmax = max(wipunit);
            pp = pmax - pmin;
            wipunit2 = wipunit./pp;
            [punit2 funit2]=periodogram(wipunit2,blackmanharris(length(wipunit2)),65536,10000);
            relthetapower_norm = (sum(punit2(21:41))) / (sum(punit2(4:end)));    %normalized relative theta power
            maxamp_norm = max(punit2(21:41));    %normalized maximum theta power
            if relthetapower > 0.25 & length(vdisc) > 10,       %1st & 2nd criterium
                %disp(['No theta modulated activity.']);
                            
% Return plot
                ivs=diff(vdisc);
                livs=length(ivs);
                d1=ivs(1:livs-1);
                d2=ivs(2:livs);
                % figure
                % plot(d1,d2,'.')
                % title('return plot of ISIs')
                % xlabel('nth ISI')
                % ylabel('(n+1)th ISI')
                % k=waitforbuttonpress;
                % if k==0 close; end;
                
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
                
                switch maxa3(1)
                case 1 
                    ivss3=ivs(c31);
                case 2
                    ivss3=ivs(c32); 
                end;
                binum1=fix(0.606+0.4*(length(ivss3)-1));
                [n1 m1]=hist(ivss3,binum1);
                K3 = find(ivss3==min(ivs));
                L3 = find(ivss3>1500);
                
                switch maxa4(1)
                case 1 
                    ivss4=ivs(c41); 
                case 2 
                    ivss4=ivs(c42);
                case 3 
                    ivss4=ivs(c43); 
                end;
                binum2=fix(0.606+0.4*(length(ivss4)-1));
                [n2 m2]=hist(ivss4,binum2);
                K4 = find(ivss4==min(ivs));
                L4 = find(ivss4>1500);
                
                switch maxa5(1)
                case 1 
                    ivss5=ivs(c51);
                case 2 
                    ivss5=ivs(c52); 
                case 3 
                    ivss5=ivs(c53); 
                case 4 
                    ivss5=ivs(c54);
                end;
                binum3=fix(0.606+0.4*(length(ivss5)-1));
                [n3 m3]=hist(ivss5,binum3);
                K5 = find(ivss5==min(ivs));
                L5 = find(ivss5>1500);
                
                switch maxa6(1)
                case 1 
                    ivss6=ivs(c61);
                case 2
                    ivss6=ivs(c62); 
                case 3 
                    ivss6=ivs(c63);
                case 4 
                    ivss6=ivs(c64);
                case 5
                    ivss6=ivs(c65);
                end;
                binum4=fix(0.606+0.4*(length(ivss6)-1));
                [n4 m4]=hist(ivss6,binum4);
                K6 = find(ivss6==min(ivs));
                L6 = find(ivss6>1500);
                
                switch maxa7(1)
                case 1
                    ivss7=ivs(c71);
                case 2
                    ivss7=ivs(c72);
                case 3
                    ivss7=ivs(c73);
                case 4 
                    ivss7=ivs(c74);
                case 5 
                    ivss7=ivs(c75);
                case 6 
                    ivss7=ivs(c76); 
                end;
                binum5=fix(0.606+0.4*(length(ivss7)-1));
                [n5 m5]=hist(ivss7,binum5);
                K7 = find(ivss7==min(ivs));
                L7 = find(ivss7>1500);
                
                [ng mg]=hist(ivs,binum);
                lmg=fix(max(mg));
                
                % figure
                % subplot(2,1,1)
                % bar(mg,ng)
                % subplot(2,1,2)
                % bar(m1,n1)
                % set(gca,'XLim',[0 lmg])
                % figure
                % subplot(2,1,1)
                % bar(mg,ng)
                % subplot(2,1,2)
                % bar(m2,n2)
                % set(gca,'XLim',[0 lmg])
                % figure
                % subplot(2,1,1)
                % bar(mg,ng)
                % subplot(2,1,2)
                % bar(m3,n3)
                % set(gca,'XLim',[0 lmg])
                % figure
                % subplot(2,1,1)
                % bar(mg,ng)
                % subplot(2,1,2)
                % bar(m4,n4)
                % set(gca,'XLim',[0 lmg])
                % figure
                % subplot(2,1,1)
                % bar(mg,ng)
                % subplot(2,1,2)
                % bar(m5,n5)
                % set(gca,'XLim',[0 lmg])
                %dec=input('Choose the No of clusters you like best -or 0 for none:-)');
                
                q = zeros(1,7);     %3rd & 4th criterium
                for i = 3:7,
                    w = [eval(['isempty(K',num2str(i),')'])... 
                            eval(['isempty(L',num2str(i),')'])...
                            eval(['isempty(L',num2str(i),')'])];
                    q(i) = sum(w);
                end;
                e = find(q==2);
                if isempty(e),
                    dec = 0;
                else dec = min(e);
                end;
                
                if dec ~= 0,
                    eval(['liv = find(ivs<max(ivss',num2str(dec),'));']);
                    liv1=[-1 liv]; liv=[liv 0];
                    bliv=find(liv~=liv1+1);
                    blivv=bliv(2:end)-1; bliv(end)=[];
                    burst=[liv(bliv);liv(blivv)+1];     %1st and last spikes of bursts
                    fsp=vdisc(burst(1,:));               %localisation of 1st spikes of bursts
                    diffburst = vdisc(burst(2,:)) - vdisc(burst(1,:));
                    mdb = max(diffburst);
                    
%                     evif = 1;
%                     while isempty(find(dec)) == 0 & evif == 1,
%                         if mdb > 2000,       %5th criterium
%                             e(dec) = 0;
%                             dec = min(e);
%                             evif = 1;
%                         else evif = 0;
%                         end;
%                     end;
                    while mdb > 3000,
                        fnd = find(diffburst == mdb);
                        burst = [burst(1,1:fnd(1)-1) burst(1,fnd(1)+1:end);...
                                burst(2,1:fnd(1)-1) burst(2,fnd(1)+1:end)];
                        diffburst = vdisc(burst(2,:)) - vdisc(burst(1,:));
                        mdb = max(diffburst);
                    end;
                end;
                
                switch dec
                case 0
                    %disp('Your cell is not bursty -such a shame:-(')
                    bursting=0;
                case 3
                     bursting=1;
                    blf=1/(max(ivss3)/10000);
                    aveintra=1/(mean(ivss3)/10000);
                    stdintra=1/(std(ivss3)/100000);
                case 4
                    bursting=1;
                    blf=1/(max(ivss4)/10000);
                    aveintra=1/(mean(ivss4)/10000);
                    stdintra=1/(std(ivss4)/100000);
                case 5
                    bursting=1;
                    blf=1/(max(ivss5)/10000);
                    aveintra=1/(mean(ivss5)/10000);
                    stdintra=1/(std(ivss5)/100000);
                case 6
                    bursting=1;
                    blf=1/(max(ivss6)/10000);
                    aveintra=1/(mean(ivss6)/10000);
                    stdintra=1/(std(ivss6)/100000);
                case 7
                    bursting=1;
                    blf=1/(max(ivss7)/10000);
                    aveintra=1/(mean(ivss7)/10000);
                    stdintra=1/(std(ivss7)/100000);
                end;
                
                burstout=struct('blf',{0},'ibsn',{0},'ibfr',{0},'blen',{0}, ...
                    'aibsn',{0},'ablen',{0}, ...
                    'sibfr',{0},'sibsn',{0},'sblen',{0});
                
                %marking of the bursts & burst characteristics
                if bursting,                    
                    bn=size(burst,2);
                    burstout.blf=blf;
                    burstout.ibsn=burst(2,:)+1-burst(1,:);
                    burstout.blen=vdisc(burst(2,:))-vdisc(burst(1,:));
                    burstout.ibfr=aveintra;
                    burstiness=sum(burstout.ibsn)/length(vdisc);
                    burstout.ablen=sum(burstout.blen)/(10*bn);
                    burstout.aibsn=sum(burstout.ibsn)/bn;
                    burstout.sblen=std(burstout.blen)/(10*bn);
                    burstout.sibfr=stdintra;
                    burstout.sibsn=std(burstout.ibsn)/bn;
                    bfs=burst(1,:);
                end;
                
                dt=0.0001;
                time=[0:length(unit)-1]*dt;
                close all
                h = figure;
%                 subplot(2,1,1);
                plot(time,unit,'m');
                if bursting == 1,
                    hold on;
                    mi = kuszob; 
                    ma = kuszob+0.5;
                    LE = zeros(1,length(burst));
                    for j = 1:length(burst),
                        rajz = [time(vdisc(burst(1,j))) time(vdisc(burst(2,j)));mi ma];
                        line(rajz(1,:),rajz(2,:),'color','b');
                        b = vdisc(burst(1,j):burst(2,j))./10;
                        eval(['B',num2str(j),'=diff(b);']);
                        eval(['LE(j)=length(B',num2str(j),');']);
                    end;
                    intraburstiv = zeros(1,sum(LE));
                    next = 1;
                    for j = 1:length(burst),
                        for k = 1:LE(j),
                            eval(['intraburstiv(next) = B',num2str(j),'(k);']);
                            next = next + 1;
                        end;
                    end;
                    intraburstivvar = var(intraburstiv);    %variance of the intraburst interval length
                    hold off;
                    
%                     sn1 = zeros(length(fsp),(datinx2-datinx1));
%                     for u = 1:length(fsp),
%                         sn1(u,1:(datinx2-datinx1-1)*10+1) = sinc(0.02*[(-1)*fsp(u)+1:0.1:(datinx2-datinx1-fsp(u))]);
%                     end;
%                     sn2 = sum(sn1);
%                     sn3 = fft(sn2,65536);
%                     sn4 = sn3 .* conj(sn3) / 512;
%                     plot(sn4(1:100));
                    
%                     iipunit = zeros(1,length(unit));
%                     iipunit(fsp) = 1;
%                     wwbh = blackmanharris(1000);
%                     wwipunit = conv(ipunit,wbh);
%                     wwipunit = wwipunit(1:length(unit));
%                     [ppunit ffunit] = periodogram(wwipunit,blackmanharris(length(wwipunit)),65536,10000);
%                     subplot(2,1,2);
%                     plot(ffunit(1:100),ppunit(1:100));
                    
                end;
                if bursting,
                    title('Bursts');
                    ymin = min(unit)-0.5;
                    ymax = max(unit)+1;
                    xlim = xlim;
                    axis([xlim(1) xlim(2) ymin ymax]);
                    ylim = ylim;
                    text(xlim(1)+1,ylim(2)-0.1,['{\itRelative theta power: }',num2str(relthetapower),'(',num2str(relthetapower_norm),')']);
                    text(xlim(1)+1,ylim(2)-0.3,['{\itMaximum theta power: }',num2str(maxamp),'(',num2str(maxamp_norm),')']);
                    text(xlim(1)+1,ylim(2)-0.5,['{\itThe frequency where theta power is maximal: }',num2str(maxfreq)]);
                    text(xlim(1)+1,ylim(2)-0.7,['{\itVariance of the intraburst interval length }',num2str(intraburstivvar)]);
                    clear ylim xlim;
                else title('No bursts');
                end;
                
% Saving the plot
                pont = findstr(filename,'.');
                filenam = filename(1:pont(1)-1);
                switch p
                case 1
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6)]);
                    eval(['saveas(h,''NO_THETA_',filenam,'.fig'')']);
                    ffilenam = [filenam,' NO THETA'];
                    fprf = [ffilenam,' ', num2str(relthetapower),' ',num2str(relthetapower_norm),' ', num2str(maxamp),...
                            ' ',num2str(maxamp_norm),' ',num2str(maxfreq),' ',num2str(intraburstivvar)];
                    fprintf(fid,'%s %s %s %s %s %s %s\n',fprf);
                    fprintf(fid,'\n',fprf);
                case 2
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6)]);
                    eval(['saveas(h,''SP_THETA_',filenam,'.fig'')']);
                    ffilenam = [filenam,' SP THETA'];
                    fprf = [ffilenam,' ', num2str(relthetapower),' ',num2str(relthetapower_norm),' ', num2str(maxamp),...
                            ' ',num2str(maxamp_norm),' ',num2str(maxfreq),' ',num2str(intraburstivvar)];
                    fprintf(fid,'%s %s %s %s %s %s %s\n',fprf);
                    fprintf(fid,'\n',fprf);
                case 3
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6)]);
                    eval(['saveas(h,''EV_THETA_',filenam,'.fig'')']);
                    ffilenam = [filenam,' EV THETA'];
                    fprf = [ffilenam,' ', num2str(relthetapower),' ',num2str(relthetapower_norm),' ', num2str(maxamp),...
                            ' ',num2str(maxamp_norm),' ',num2str(maxfreq),' ',num2str(intraburstivvar)];
                    fprintf(fid,'%s %s %s %s %s %s %s\n',fprf);
                    fprintf(fid,'\n',fprf);   
                end;
            end;
        end;
    end;
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
fclose(fid);
cd(mm);