%% import

loadcb
fullpth = 'c:\Balazs\_anatomy\NB\alltracks_new.xlsx';
[t0 t1 cdata0] = xlsread(fullpth,'session_table');     % Excel
cdata = cdata0;
for k = 1:size(cdata0,1)
    cdata{k,2} = cdata{k,2}(1:end-1);
end

%% match

NumCells = length(CELLIDLIST);
newdata = cell(NumCells,6);
for k = 1:NumCells
    cellid = CELLIDLIST{k};
    [a s] = cellid2tags(cellid);
    if isequal(a,'n007')    % no track reconstruction
        continue
    end
    sinx = find(cellfun(@(v)~isempty(v),strfind(cdata(:,2),s)));
    ainx = find(cellfun(@(v)~isempty(v),strfind(cdata(:,1),a)));
    inx = intersect(ainx,sinx);
    if length(inx) ~= 1
        error('No unique match.')
    end
    newdata{k,1} = cellid;
    newdata(k,2:6) = cdata(inx,3:7);
end

%% import - n023 TT1

loadcb
fullpth = 'c:\Balazs\_anatomy\NB\alltracks_new.xlsx';
[t0 t1 cdata0_n023] = xlsread(fullpth,'n023_TT1');     % Excel
cdata_n023 = cdata0_n023;
for k = 1:size(cdata0_n023,1)
    cdata_n023{k,2} = cdata_n023{k,2}(1:end-1);
end

%% match - special case for n023 TT1

NumCells = length(CELLIDLIST);
for k = 1:NumCells
    cellid = CELLIDLIST{k};
    [a s t] = cellid2tags(cellid);
    if ~isequal(a,'n023') || ~isequal(t,1)
        continue
    end
    sinx = find(cellfun(@(v)~isempty(v),strfind(cdata_n023(:,2),s)));
    ainx = find(cellfun(@(v)~isempty(v),strfind(cdata_n023(:,1),a)));
    inx = intersect(ainx,sinx);
    if length(inx) ~= 1
        error('No unique match.')
    end
    newdata{k,1} = cellid;
    newdata(k,2:6) = cdata_n023(inx,3:7);
end

%% import - n028

loadcb
fullpth = 'c:\Balazs\_anatomy\NB\alltracks.xlsx';
[t0 t1 cdata0_n028_2_7] = xlsread(fullpth,'n028_TT2_TT7');     % Excel
cdata_n028_2_7 = cdata0_n028_2_7;
for k = 1:size(cdata0_n028_2_7,1)
    cdata_n028_2_7{k,2} = cdata_n028_2_7{k,2}(1:end-1);
end

[t0 t1 cdata0_n028_others] = xlsread(fullpth,'n028_others');     % Excel
cdata_n028_others = cdata0_n028_others;
for k = 1:size(cdata0_n028_others,1)
    cdata_n028_others{k,2} = cdata_n028_others{k,2}(1:end-1);
end

%% match - special case for n028

NumCells = length(CELLIDLIST);
for k = 1:NumCells
    cellid = CELLIDLIST{k};
    [a s t] = cellid2tags(cellid);
    if isequal(a,'n028') && ismember(t,[2 7])
        sinx = find(cellfun(@(v)~isempty(v),strfind(cdata_n028_2_7(:,2),s)));
        ainx = find(cellfun(@(v)~isempty(v),strfind(cdata_n028_2_7(:,1),a)));
        inx = intersect(ainx,sinx);
        if length(inx) ~= 1
            error('No unique match.')
        end
        newdata{k,1} = cellid;
        newdata(k,2:6) = cdata_n028_2_7(inx,3:7);
    elseif isequal(a,'n028') && ~ismember(t,[2 7])
        sinx = find(cellfun(@(v)~isempty(v),strfind(cdata_n028_others(:,2),s)));
        ainx = find(cellfun(@(v)~isempty(v),strfind(cdata_n028_others(:,1),a)));
        inx = intersect(ainx,sinx);
        if length(inx) ~= 1
            error('No unique match.')
        end
        newdata{k,1} = cellid;
        newdata(k,2:6) = cdata_n028_others(inx,3:7);
    end
end

%% A1

NumCells = length(CELLIDLIST);
for k = 1:NumCells
    cellid = CELLIDLIST{k};
    [a s t] = cellid2tags(cellid);
    if (isequal(a,'n013') && isequal(t,6)) ||...
            (isequal(a,'n018') && isequal(t,5)) ||...
            (isequal(a,'n020') && isequal(t,8)) ||...
            (isequal(a,'n021') && isequal(t,5)) ||...
            (isequal(a,'n023') && isequal(t,5)) ||...
            (isequal(a,'n026') && isequal(t,5)) ||...
            (isequal(a,'n027') && isequal(t,5)) ||...
            (isequal(a,'n028') && isequal(t,5)) ||...
            (isequal(a,'n029') && isequal(t,5)) ||...
            (isequal(a,'n037') && isequal(t,8)) ||...
            (isequal(a,'n038') && isequal(t,8)) ||...
            (isequal(a,'n039') && isequal(t,8)) ||...
            (isequal(a,'n040') && isequal(t,5)) ||...
            (isequal(a,'n043') && isequal(t,5)) ||...
            (isequal(a,'n045') && isequal(t,5)) ||...
            (isequal(a,'n046') && isequal(t,5))
        
        newdata{k,1} = cellid;
        newdata(k,2:6) = {2400 -2500 4000 'ACx' 'ACx'};
    end
end

%% save

xlswrite(fullpth,newdata,'cell_table','A1')