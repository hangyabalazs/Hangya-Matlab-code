% Imported from Viktor.

% Same as multi, just skip the cell where soome problem was during the analysis
% Needs multicatchfcn: (nf=nf+1); 
%       multi_f:the analysis part

dont='0';
if exist('inputpath') & exist('figmentpath') & exist('paramentpath'),
    disp(['Input path is:             ' inputpath]);
    disp(['Figure saving path is:     ' figmentpath]);
    disp(['Parameters saving path is: ' paramentpath]);    
    dont=input('Are theese correct? /ENTER: yes; 0: no/ ');
end;
if ~isempty(dont),
[fname inputpath]=uigetfile( ...
    {'*.*',  'All Files (*.*)'; ...
     '*.txt', 'ASCII Files (*.txt)'; ...
     '*.mat', 'Mat Files (*.mat)'}, ...
     'Sign a file in that folder where the raw data files are!');

[fname figmentpath]=uigetfile( ...
    {'*.*',  'All Files (*.*)'; ...
     '*.txt', 'ASCII Files (*.txt)'; ...
     '*.mat', 'Mat Files (*.mat)'}, ...
     'Sign a file in that folder where figures will be saved!');
[fname paramentpath]=uigetfile( ...
    {'*.*',  'All Files (*.*)'; ...
     '*.txt', 'ASCII Files (*.txt)'; ...
     '*.mat', 'Mat Files (*.mat)'}, ...
     'Sign a file in that folder where the cell parameters will be saved!');
end;

fnames=dir(inputpath);

nf=1; notanal=''; kuszobtpe=0; kuszobtpu=0;
while nf<size(fnames,1)+1,
    disp(['Inporting data ...']);
    clear fname;
    while ~exist('fname') & nf<size(fnames,1)+1,
        if ~fnames(nf).isdir, 
            fname=[inputpath fnames(nf).name]; 
            
            fn=[paramentpath '\' fnames(nf).name(1:6) '.txt'];
            if ~exist(fn),
                fpara=fopen(fn,'w');
            else 
                disp(['There is a .txt file for ' fnames(nf).name(1:6) ' data at the output directory! ']);
                disp(' Do you want to analyse the cell or not?');
                dont=input(' /Yes: Enter; No: 0/ ');
                if ~isempty(dont),
                    notanal=[notanal fnames(nf).name(1:6) ' '];
                    clear fname;
                    nf=nf+1;
                else 
                    disp('Are the saved parameters acceptables? ')
                    dont=input(' /Yes: Enter; No: 0/ ');
                    if isempty(dont),
                        clear fname;
                        nf=nf+1;
                    else 
                        fpara=fopen(fn,'w');
                    end;
                end;
            end;
        else 
            nf=nf+1;
        end;
    end;
                    
  if nf<size(fnames,1)+1,
    disp(['Inporting data ' fnames(nf).name(1:6) ' ...']);
    data=load(fname);
        
    if isstruct(data)
        field=fieldnames(data);
        s=struct('type','.','subs',field);
        data=subsref(data,s);
    end;
    fname=fnames(nf).name;
    disp([fname(1:6) ' is inported.']);
    
    if rem(length(data),2)==1
        data(length(data)+1)=0;
    end;
    if size(data,2)==1, data=[data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)]; end;
    
    disp('Press ''0'' if there is no such type of activity');
    inx=zeros(1,6);
    inx(1)=input('First point of interval OF NO THETA         = ');
    inx(2)=input('Last point of interval OF NO THETA          = ');
    inx(3)=input('First point of interval OF SPONTAN THETA    = ');
    inx(4)=input('Last point of interval OF SPONTAN THETA     = ');
    inx(5)=input('First point of interval OF EVOKED THETA     = ');
    inx(6)=input('Last point of interval OF EVOKED THETA      = ');
    
    fprintf(fpara,'inxnoth1    %11e\n',inx(1));
	fprintf(fpara,'inxnoth2    %11e\n',inx(2));
	fprintf(fpara,'inxspth1    %11e\n',inx(3));
	fprintf(fpara,'inxspth2    %11e\n',inx(4));
	fprintf(fpara,'inxevth1    %11e\n',inx(5));
    fprintf(fpara,'inxevth2    %11e\n',inx(6));
    
    i=1;
    while i<length(inx)+1,
        if inx(i)==0, inx(i)=[]; i=i-1; end;
        i=i+1;
    end;
    
    for i=1:2:length(inx),
	    datinx1=inx(i); datinx2=inx(i+1);
        unit=data(datinx1:datinx2,2)'; 
        eeg=data(datinx1:datinx2,1)'; 
        dt=0.0001;
        time=[0:length(unit)-1].*dt;
       
        s=figure;
        plot(time,unit,'m');
        set(gcf,'numbertitle','off', ...
                'name',[fname(1:6) '_' num2str(datinx1/10000) '_' num2str(datinx2/10000)]);
        hold on
        kuszob=2*(max(unit)-min(unit))/3+min(unit);
        plot([time(1) time(end)],[kuszob kuszob]);
        title('The treshold remains: ENTER /Command window/')
        disp('The treshold?');
        disp('/ENTER: the treshold remains; ANY OTHER NUMBER: the new treshold/ ');
        ok=0;
        while ~ok, 
            dont=input('');
            if isempty(dont), 
                disp(['The treshold is ' num2str(kuszob)]);
                ok=1;
            elseif  isnumeric(dont), 
                kuszob=dont;
                ok=1;
            else 
                disp(['You must tape a number! Tape it or ENTER for treshold remains ' num2str(kuszob)]);
            end;
        end;
        hold off
        fprintf(fpara,'kuszob    %11e\n',kuszob);
        close(s);
    end;
	fclose(fpara); 
    nf=nf+1;
  end;
end;

nf=1;
fnames=dir(inputpath);
lasterr('');
while nf<size(fnames,1)+1;
sok_f2b
    %     eval('sok_f2b','multicatchfcn');
end;