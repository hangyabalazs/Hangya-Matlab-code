function b_empty2004(RunFile,DataDir,ResultDir,IntervalSet,TextOutput,Name,InputArg)
%EMPTY2004   Batch processing on hippocampus data files.
%   Different from EMPTY and EMPTY_NEW, you have to specify the parameters via EMPTY_GUI2004
%   Graphical User Interface. Start EMPTY_GUI2004 to run!
%
%   See also EMPTY, EMPTY_NEW and EMPTY_GUI2004.

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;

where1 = DataDir;    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd(ResultDir);  %Here are the results
if ~isempty(TextOutput)
    npt = input(['Discard existing content or append data while writing ',...
            TextOutput,'.txt? /discard:  ENTER, append: a/'],'s');
    if isempty(npt)
        str = ['fid = fopen(' TextOutput ',''w'');'];
        eval(str);
    elseif npt == 'a'
        str = ['fid = fopen(' TextOutput ',''a'');'];
        eval(str);
    else
        error('Unexpected answer for input.')
    end
    global FID
    FID = fid;
end

% Interval set
if isequal(IntervalSet,'manual')
    datinx1 = input('Datinx1: ');
    datinx2 = input('Datinx2: ');
end

% Load data
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end
    if isequal(IntervalSet,'auto')
        datinx1 = 1;
        datinx2 = size(data,1);
    end
    meret = size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end
    if size(data,2) == 1
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
            
%Computing the input variables
    dt = 0.0001;
    mintafr = 1 / dt;
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
    IN{4} = where1;
    IN{5} = datinx1;
    IN{6} = datinx2;
    IN{7} = time;
    IN{8} = unit;
    IN{9} = dt;
    IN{10} = meret;
    IN{11} = mintafr;
    IN{12} = xlimit;
            
% % Discrimination                
%     b_disc(kuszob);
%     global DISC
%     id = DISC{1};
%     output = DISC{2};
%     vdisc = DISC{3};
%     kuszob = DISC{4};
%     instfrek = DISC{5};
%     isi = DISC{6};
            
% Main
            str = ['OutputArg = ',RunFile,'(InputArg);'];
            eval(str);

% Saving
            if ~isstruct(OutputArg)
                error('Output argument of the main function has to be structure array.')
            end
            vgo = fieldnames(OutputArg);
            for s = 1:length(vgo);
                vg = vgo{s};
                str = ['v = OutputArg.' vg ';'];
                eval(str);
                if ishandle(v)
                    eval(['saveas(v,''',Name,'_',fname(1:6),'.fig'')']);
                else
                    eval([vg ' = v;'])
                    eval(['save(''',Name,'_',fname(1:6),'.mat'',''',vg,''')']);
                end
            end
        end
    end
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator

if ~isempty(TextOutput)
    fclose(fid);
end

close all
cd(mm);