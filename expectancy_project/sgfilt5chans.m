% filename={'norbi_corr.set','hg06corr.set','hg07corr.set','hg08corr.set',...
%     'hg09corr.set','hg10corr.set','hg11corr.set', 'hg1201.set_concat.set',...
%     'hg1301.set_concat.set', 'hg1401.set_concat.set', 'hg1501.set_concat.set'...
%     'hg1601.set_concat.set', 'hg1701.set_concat.set'}; cd('D:\ATTENTION_HIGHGAMMA\setdataconcat');
filename=dir('X:\In_Vivo\raw_data\human_SG\setdataconcatFZCZPZC3C4/*.set'); 
filename=struct2cell(filename); 
filename=filename(1,:);
for f=1:size(filename,2)

EEG=pop_loadset(strcat('X:\In_Vivo\raw_data\human_SG\setdataconcatFZCZPZC3C4/', char(filename(f)))); 
% cd('D:\ATTENTION_HIGHGAMMA\setdataconcat'); D:\ATTENTION_HIGHGAMMA\setdataconcatFZCZPZ
%EEG = pop_select( EEG, 'channel',[2 4 5 6 8] ); %EEG = pop_select( EEG, 'channel',[2:3:8] ); 
EEG = eeg_checkset( EEG ); % 2:3:

%%FIR1 filtering (matlab built-i)
fprintf('..FIR1 Filtering data...');
flt = fir1(2048,[0.5 3]/500); %4096
for m = 1:size(EEG.data,1)
EEG.data(m,:) = filtfilt(flt,1,EEG.data(m,:));
end
fprintf('...done...');

save(strcat('D:\ATTENTION_HIGHGAMMA\setdataconcatFZCZPZC3C4filt/', char(filename(f))), 'EEG');
end
