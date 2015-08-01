%% load attention cells

load('C:\Balazs\_analysis\NB\attentioncells3\summary\attentioncells.mat')

%% reconstruct attention cells

reconstruction(cellids)

%% FR

ATTbaseline = getvalue('Baseline',cellids);

%% Spike shape

ATTspikeshape = getvalue('SpikeShape',cellids);
ATTspikeshape = nancell2struct(ATTspikeshape);
ATTspike = {ATTspikeshape.Spike};
ATTspike = cell2mat(cellfun(@(s)zscore(s'),ATTspike,'UniformOutput',false))';

%% load attention cells

load('C:\Balazs\_analysis\NB\accuracycells\summary\accuracycells.mat')

%% reconstruct attention cells

reconstruction(cellids)

%% FR

ACCbaseline = getvalue('Baseline',cellids);

%% Spike shape

ACCspikeshape = getvalue('SpikeShape',cellids);
ACCspikeshape = nancell2struct(ACCspikeshape);
ACCspike = {ACCspikeshape.Spike};
ACCspike = cell2mat(cellfun(@(s)zscore(s'),ACCspike,'UniformOutput',false))';