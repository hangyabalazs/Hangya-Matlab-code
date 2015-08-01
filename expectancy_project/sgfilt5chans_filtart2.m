%SGFILT5CHANS_FILTART2   Filtering (from S.G.).
%   SGFILT5CHANS_FILTART2 creates filtered data from raw data of
%   Experiment2. It requires the EEGlab Matlab package. Subsequent programs
%   are SGEPOCH2, SGEPOCHPHASE2 and SGCONTROLPHASE2.
%
%   See also SGEPOCH2, SGEPOCHPHASE2 and SGCONTROLPHASE2.

% filename={'norbi_corr.set','hg06corr.set','hg07corr.set','hg08corr.set',...
%     'hg09corr.set','hg10corr.set','hg11corr.set', 'hg1201.set_concat.set',...
%     'hg1301.set_concat.set', 'hg1401.set_concat.set', 'hg1501.set_concat.set'...
%     'hg1601.set_concat.set', 'hg1701.set_concat.set'}; cd('D:\ATTENTION_HIGHGAMMA\setdataconcat');
filename=dir('X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2/*.set'); 
filename=struct2cell(filename); 
filename=filename(1,:);
for f=1:size(filename,2)

EEG=pop_loadset(strcat('X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2/', char(filename(f)))); 
% cd('D:\ATTENTION_HIGHGAMMA\setdataconcat'); D:\ATTENTION_HIGHGAMMA\setdataconcatFZCZPZ
%EEG = pop_select( EEG, 'channel',[2 4 5 6 8] ); %EEG = pop_select( EEG, 'channel',[2:3:8] ); 
EEG = eeg_checkset( EEG ); % 2:3:

%%FIR1 filtering (matlab built-i)
fprintf('..FIR1 Filtering data...');
flt = fir1(2048,[0.5 3]/250); %4096
for m = 1:size(EEG.data,1)
% m = 3;
eeg1 = EEG.data(m,:);
for k = 1:length(EEG.event)
    eet = EEG.event(k).type;
    if isequal(eet,'22') || isequal(eet,'42')
        eel = round(EEG.event(k).latency);
        einx = min(eel+500,length(eeg1));
        einxx = eel+1:einx;
        eeg1(einxx) = randn(size(einxx)) + eeg1(eel);
    end
end
EEG.data(m,:) = filtfilt(flt,1,eeg1);
% eeg1 = filtfilt(flt,1,eeg1);
end
fprintf('...done...');

EEG_Cz = eeg1;
save(strcat('X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2filtart/', char(filename(f))), 'EEG');
end
