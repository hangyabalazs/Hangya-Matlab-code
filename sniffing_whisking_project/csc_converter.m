function csc_converter(fn)

% Convert CSC files
for k = 1:8
    csc = ['CSC' num2str(k) '.ncs'];
    [TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
        Nlx2MatCSC([fn csc],[1 1 1 1 1],1,1,1);  % convert CSC
    fns = fullfile(fn,['CSC' num2str(k) '.mat']);
    save(fns,'TimeStamp','ChanNum','SampleFrequency','NumValSamples','Samples','NlxHeader')
    
end

% Convert Event file
[EventTimeStamps, EventIDs, Nttls, Extras, EventStrings NlxHeader] = ...
    Nlx2MatEV([fn 'Events.nev'],[1 1 1 1 1],1,1,1);
fns2 = fullfile(fn,'Events.mat');
save(fns2,'EventTimeStamps','EventIDs','Nttls','Extras','EventStrings','NlxHeader');