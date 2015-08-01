function ichanmatch

trcfile = 'F:\raw_data\human_SO\oiti37_lukacs\TRC\patients\pat_1\EEG_12.trc';
matfile = 'F:\raw_data\human_SO\oiti37_lukacs\grid\mat\EEG12_4x5grid_18i.mat';
load(matfile)
fid = fopen(trcfile);
trcdata = trc_import(fid);
sd = size(data,2);
chan = zeros(1,sd);
for k = 1:sd
    datarow = data(:,k);
    trm = 0;
    next = 1;
    while ~trm
        trcrow = trcdata;
        if trcrow - datarow < 0.1
            chan(k) = next;
            trm = 1;
        end
        next = next + 1;
    end
end

% -------------------------------------------------------------------------
function trcdata = trc_import(f)

fseek(f,175,-1);            % headertype
d = fread(f,1,'int8');
switch d
    case 1
        fseek(f,182,-1);
        minbyte = fread(f,1,'int16');
        orig_chnum = fread(f,1,'int16');
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
        end
    case 4
        fseek(f,138,-1);
        minbyte = fread(f,1,'int32');
        orig_chnum = fread(f,1,'int16');
        multiplexer = fread(f,1,'int16');
        srate = fread(f,1,'int16');
        databyte = fread(f,1,'int16');
        if databyte * orig_chnum ~= multiplexer,
            msgbox('Multiplexer error');
        end
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
        end
    case {2,3}
        fseek(f,138,-1);
        minbyte = fread(f,1,'int32');
        orig_chnum = fread(f,1,'int16');
        multiplexer = fread(f,1,'int16');
        srate = fread(f,1,'int16');
        databyte = fread(f,1,'int16');
        if databyte * orig_chnum ~= multiplexer,
            msgbox('Multiplexer error');
        end
        fseek(f,184,-1);
        ordpoi = fread(f,1,'int32');
        fseek(f,200,-1);
        elpoi = fread(f,1,'int32');
        fseek(f,ordpoi,-1);
        order = fread(f,orig_chnum,'int8');
        for a = 1:orig_chnum,
            fseek(f,elpoi+2+order(a)*128,-1);
            p = fread(f,6,'int8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,1} = cp;
            p = fread(f,6,'int8'); 
            p(p==0) = 32;
            cp = char(p');
            orig_chnames{a,2} = cp;
        end
    case 0
        disp('Type 0 header not implemented yet');
        return
end