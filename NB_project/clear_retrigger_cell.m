% load('C:\Balazs\_data\NB\NB_cellbase\n018\111014a\TT2_1.mat')
% Nttfile = 'c:\Balazs\_data\NB\NB_cellbase\n018\111014a\TT2.ntt';
% [junkTS,waveform] = LoadTT_NeuralynxNT(Nttfile,TS,1);
% 
% sinx = find(diff(TS/10)<0.75);

% Input directory
inpdir = 'c:\Balazs\_data\NB\NB_cellbase\n018\111015a\';
dr = dir(inpdir);
for k = 3:length(dr)
    
    % Select
    fn = dr(k).name;
    if length(fn) < 3 || ~isequal(fn(end-3:end),'.mat') || ~isequal(fn(1:2),'TT')
        continue
    end
    disp(fn)
    
    % Load
    fname = fullfile(inpdir,fn);
    load(fname)
    
    % Drop dupliacate spikes
    sinx = find(diff(TS/10)<0.75);
    TS(sinx+1) = [];
    
    % Save
    save(fname,'TS')
end