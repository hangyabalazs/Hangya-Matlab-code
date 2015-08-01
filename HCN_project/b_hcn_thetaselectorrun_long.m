function b_hcn_thetaselectorrun_long
%HCN_THETASELECTORRUN_LONG   Version of HCN_THETASELECTORRUN for long registrations.
%   You have to use HCN_THETASELECTORRUN_LONG for registrations longer than 200 seconds.
%   HCN_THETASELECTORRUN saves a binary file called 'long.mat', which contains the 
%   identifiers of the long registrations. You have to put the "long files" in the 
%   input directory of HCN_THETASELECTORRUN_LONG (either manually or using LONGCOPY).
%   For details, see HCN_THETASELECTORRUN.
%
%   HCN_THETASELECTORRUN_LONG uses two directories: one for the data files and one for
%   the results. You are able to modify these directories through editing the program
%   code.
%
%   See also THETASELECTOR3, THETASELECTOR_BETA3, LONGCOPY, THETASELECTORRUN,
%   THETASELECTORRUN_LONG and HCN_THETASELECTORRUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
global DATADIR
if isempty(DATADIR)
    clear glogal DATADIR;
    b_startup
    global DATADIR
end;
% where1 = DATADIR;    %Here are the data files
where1 = ['f:\raw_data\hcn\temp\'];
files = dir(where1);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
    end
end
files2 = files2(2:end);
sf = length(files2);
mm = pwd;
cd([DATAPATH, 'HCN\Wavelet2\thetaselection_long\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    for dtx = 1:2
        fname = files2(o).name;
        ffnm = [where1 fname];
        data = load(ffnm);
        meret=size(data,1);
        if isstruct(data),
            field = fieldnames(data);
            s = struct('type','.','subs',field);
            data = subsref(data,s);
        end;
        if size(data,2)==1,
            data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
        end;
        
% Computing the input variables
        if dtx == 1
            datinx1 = 1; %first point of the interval
            datinx2 = 2000000;  %last point of the interval
        else
            datinx1 = 2000001;  %first point of the interval
            datinx2 = size(data,1); %last point of the interval
        end
        if datinx2 - datinx1 > 2000000
            global LONG
            LONG{end+1} = fname(1:6);
            continue
            waitbar(o/sf)   %Progress indicator
        end
        
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
        IN{4} = where1;     % pathname
        IN{5} = datinx1;
        IN{6} = datinx2;
        IN{7} = time;
        IN{8} = unit;
        IN{9} = dt;
        IN{10} = meret;
        IN{11} = mintafr;
        IN{12} = xlimit;
        
% Free memory
        save temp o sf wb files2 where1 dtx
        clear
        clear function
        
% Theta selection
        [wave,f,newstep] = b_waveletcall_for_applications;
        save newstep newstep
        save scalefreq f
        
        global IN       % free memory
        fname = IN{3};
        clear global IN
        
        sw1 = size(wave,1);
        sw2 = size(wave,2);
        pieceno = 10;
        segm = fix(sw2/pieceno);
        inx2 = 0;
        next = 1;
        while inx2 < sw2
            inx1 = max(1,inx2+1);
            inx2 = min(sw2,inx1+segm);
            wavefrag = wave(:,inx1:inx2);
            str = ['save temp' num2str(next) ' wavefrag'];
            eval(str)
            clear wavefrag
            next = next + 1;
        end
        clear wave
        power = [];
        phase = [];
        for wsn = 1:next-1
            str = ['load temp' num2str(wsn)];
            eval(str)
            powerfrag = (abs(wavefrag)) .^2;
            power = [power powerfrag];
            clear powerfrag
            phasefrag = angle(wavefrag);
            phase = [phase phasefrag];
            clear phasefrag
            clear wavefrag
        end        
        
%         sw2 = size(wave,2);
%         pieceno = 20;
%         segm = fix(sw2/pieceno);
%         power = [];
%         while ~isempty(wave)
%             index1 = 1;
%             index2 = min(segm,size(wave,2));
%             wavefrag = wave(:,index1:index2);
%             powerfrag = (abs(wavefrag)) .^ 2;
%             clear wavefrag
%             wave(:,index1:index2) = [];
%             power = [power powerfrag];
%         end
        clear wave
        
        [H,OM] = b_thetaselector3(power,f,newstep);
        ratio = b_thetaselector_beta3(power,f,newstep);
        
% Saving
        Out = zeros(3,sw2);
        Out(1:2,:) = OM;
        Out(3,:) = ratio;
        
        pont = findstr(fname,'.');
        filenam = fname(1:pont(1)-1);
        fig = gcf;
        ax = findobj(fig,'Type','axes');
        axes(ax(2))
        title(['EEG ',filenam(1:3),' ',filenam(5:6)]);
        load temp
        if dtx == 1
            eval(['saveas(H,''THETA_SELECT_LONG1_',filenam(1:6),'.jpg'')']);
            eval(['save(''THETA_SELECT_LONG1_',filenam(1:6),'.mat'',''Out'')']);
        else
            eval(['saveas(H,''THETA_SELECT_LONG2_',filenam(1:6),'.jpg'')']);
            eval(['save(''THETA_SELECT_LONG2_',filenam(1:6),'.mat'',''Out'')']);
        end
        close all
    end
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
delete temp.mat
close all

global LONG
long = LONG;
save long long
clear global LONG