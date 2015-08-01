function iconvert3
%ICONVERT3   Converts trc to mat.
%   ICCONVERT3 loads OITI trc files and saves mat files.
%
%   For OITI39_gaal.
%
%   See also ICHANMATCH.

% Load
trcfile = 'D:\human_SO\OITI_39_GRID_02\gaalM_40.cnt';
fid = fopen(trcfile);

% Segment boundaries
beg = input(['Segment start: ']);
en = input(['Segment end: ']);
len = en - beg;

% Convert
[trcdata srate] = cnt_import(fid,beg,len);

% Display sampling rate
disp(['Sampling rate: ' num2str(srate) ' Hz'])

% Rearrange channels
K = [18:21 26:29 34:37 42:45];
KK = [46:49 54:57 66:69 74:77];
data = zeros(size(trcdata,1),16);
for k = 1:16
    data(:,k) = trcdata(:,KK(k));
end

% Save
fn = [trcfile(1:end-4) '_' num2str(beg) '_' num2str(en) '.mat'];
save(fn,'data','srate')

% -------------------------------------------------------------------------
function [data, chnames, srate] = trc_import(f,beg,len)

fseek(f,175,-1);            % headertype
d = fread(f,1,'int8');
switch d
    case 1
        fseek(f,182,-1);
        minbyte = fread(f,1,'int16');
        orig_chnum = fread(f,1,'int16');
        page = {[1:orig_chnum]};
        fseek(f,8,-1);
        srate = fread(f,1,'int16');
        databyte = 1;
        fseek(f,0,1);
        maxbyte = ftell(f);
        inx1 = 0;
        maxsec = (maxbyte - minbyte) / (orig_chnum * databyte * srate);
        fseek(f,0,-1);
        header = fread(f,minbyte,'int8');
        fseek(f,12,-1);
        mcs_factor = fread(f,1,'int16');
        fseek(f,5056,-1);
        for a = 1:orig_chnum,
            p = fread(f,5,'uint8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,1} = cp;
            orig_chnames{a,2} = 'G2';
            fread(f,1,'int8');
            d = fread(f,1,'int16');
            orig_unit(a) = fix(d/2^14);
            orig_mcs(a) = mod(d,2^14);
        end
    case 4
        fseek(f,138,-1);
        minbyte = fread(f,1,'int32');
        orig_chnum = fread(f,1,'int16');
        page = {[1:orig_chnum]};
        multiplexer = fread(f,1,'int16');
        srate = fread(f,1,'int16');
        databyte = fread(f,1,'int16');
        if databyte * orig_chnum ~= multiplexer,
            msgbox('Multiplexer error');
        end
        fseek(f,0,1);
        maxbyte = ftell(f);
        inx1 = 0;
        maxsec = (maxbyte - minbyte) / (orig_chnum * databyte * srate);
        fseek(f,0,-1);
        header = fread(f,minbyte,'int8');
        fseek(f,184,-1);
        ordpoi = fread(f,1,'int32');
        fseek(f,200,-1);
        elpoi = fread(f,1,'int32');
        fseek(f,ordpoi,-1);
        order = fread(f,orig_chnum,'int16');
        for a = 1:orig_chnum
            fseek(f,elpoi+2+order(a)*128,-1);
            p = fread(f,6,'int8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,1} = cp;
            p = fread(f,6,'int8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,2} = cp;
            orig_logic_min(a) = fread(f,1,'uint32');
            orig_logic_max(a) = fread(f,1,'uint32');
            orig_logic_gnd(a) = fread(f,1,'uint32');
            orig_physic_min(a) = fread(f,1,'int32');
            orig_physic_max(a) = fread(f,1,'int32');
            orig_unit(a) = fread(f,1,'uint16');
            orig_pref_hp(a) = fread(f,1,'uint16');
            orig_pref_hpty(a) = fread(f,1,'uint16');
            orig_pref_lp(a) = fread(f,1,'uint16');
            orig_pref_lpty(a) = fread(f,1,'uint16');
            orig_freq_coef(a) = fread(f,1,'uint16');
        end
    case {2,3}
        fseek(f,138,-1);
        minbyte = fread(f,1,'int32');
        orig_chnum = fread(f,1,'int16');
        page = {[1:orig_chnum]};
        multiplexer = fread(f,1,'int16');
        srate = fread(f,1,'int16');
        databyte = fread(f,1,'int16');
        if databyte * orig_chnum ~= multiplexer,
            msgbox('Multiplexer error');
        end
        fseek(f,0,1);
        maxbyte = ftell(f);
        inx1 = 0;
        maxsec = (maxbyte - minbyte) / (orig_chnum * databyte * srate);
        fseek(f,0,-1);
        header = fread(f,minbyte,'int8');
        fseek(f,184,-1);
        ordpoi = fread(f,1,'int32');
        fseek(f,200,-1);
        elpoi = fread(f,1,'int32');
        fseek(f,ordpoi,-1);
        order = fread(f,orig_chnum,'int8');
        for a = 1:orig_chnum
            fseek(f,elpoi+2+order(a)*128,-1);
            p = fread(f,6,'int8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,1} = cp;
            p = fread(f,6,'int8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,2} = cp;
            orig_logic_min(a) = fread(f,1,'uint16');
            orig_logic_max(a) = fread(f,1,'uint16');
            orig_logic_gnd(a) = fread(f,1,'uint16');
            orig_physic_min(a) = fread(f,1,'int32');
            orig_physic_max(a) = fread(f,1,'int32');
            orig_unit(a) = fread(f,1,'uint16');
            orig_pref_hp(a) = fread(f,1,'uint16');
            orig_pref_hpty(a) = fread(f,1,'uint16');
            orig_pref_lp(a) = fread(f,1,'uint16');
            orig_pref_lpty(a) = fread(f,1,'uint16');
            orig_freq_coef(a) = fread(f,1,'uint16');
        end
    case 0
        disp('Type 0 header not implemented yet');
        return
end

lengt = len;
begin = beg;
beginpoint = fix(begin*srate);
begin = fix(begin*1000)/1000;
lengt = fix(lengt*1000)/1000;
inlengt = lengt;
inbegin = begin;
ok = 0;
begin = round(minbyte+beginpoint*orig_chnum*databyte);
lengt = lengt * orig_chnum * srate;
fseek(f,begin,'bof');
data = fread(f,lengt,'uint16');
if orig_chnum > 1,
    n = size(data,1);
    nchn = fix(n/orig_chnum);
    data(nchn*orig_chnum+1:end) = [];
    data = reshape(data,orig_chnum,nchn)';
end
fclose(f);
chnames = orig_chnames(:,1);

% -------------------------------------------------------------------------
function [data, srate] = cnt_import(f,beg,len)

databyte = 2;
fseek(f,376,-1);
srate = fread(f,1,'int16');
fseek(f,370,-1);
chnum = fread(f,1,'int16');
page = {[1:chnum]};
if srate > 10000,
    lengt = 1;
    amp = 0.1;
else
    lengt = 10;
    amp = 0.001;
end
fseek(f,886,-1);
maxbyte = fix(fread(f,1,'int32'));
minbyte = 900 + 75 * chnum;
fseek(f,0,-1);
header = fread(f,minbyte,'int8');
inx1 = 0;
maxsec = (maxbyte - minbyte) / (chnum * databyte * srate);

lengt = len;
begin = beg;
beginpoint = fix(begin*srate);
begin = fix(begin*1000)/1000;
lengt = fix(lengt*1000)/1000;
inlengt = lengt;
inbegin = begin;
ok = 0;

begin = round(minbyte+beginpoint*chnum*databyte);
lengt = lengt * chnum * srate;
fseek(f,begin,'bof');
data = fread(f,lengt,'int16');
if chnum > 1,
    n = size(data,1);
    nchn = fix(n/chnum);
    data(nchn*chnum+1:end) = [];
    data = reshape(data,chnum,nchn)';
end;
fclose(f);
ok=1;