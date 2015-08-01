function psth = apsth2(cellid)
%APSTH2   Spike dendity function with adaptive Gaussian kernel.
%	PSTH_ACONV = APSTH2(CELLID) computes the spike density
%	function with an adaptive Gaussian kernel optimized for local
%	probability of spiking. For near-zero spiking probability, the kernel
%	is flat resulting in a moving average. For a probability of 1, the
%	kernel konverges to a Dirac delta. Note that the mapping of
%	probabilities to kernel width is arbitrary between the extremes. In
%	this implementation, the ALPHA parameter of the GAUSSWIN function
%	(which is inversly related to the SD of the Gaussian window) is
%	linearly mapped to probabilities.
%
%   See also RUN_MAKEAPSTH, BINRASTER2APSTH and GAUSSWIN.

%   Edit log: BH 8/18/2011

% Set input arguments for STIMES2BINRASTER
win = [-0.6 0.6];
dt = 0.001;   % resolution of bin raster in s
EventName = 'PulseOn';
ST = loadcb(cellid,'STIMSPIKES');    % load variables from CellBase
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName);
if epoch_pos1 == 0
    error('Epoch name not found');
end

% Time variables and valid trials
stimes = ST.event_stimes{epoch_pos1};
time = win(1):dt:win(end);
valid_trials = find(~isnan(getfield(TE,EventName)));

% Calculate bin raster
binraster = stimes2binraster(stimes(valid_trials),time,dt);

% Trial number and epoch length
[tno tl] = size(binraster);

% Calculate adaptive SDF with variable Gaussian Kernel #2
agvd = zeros(tno,tl);
prob = nanmean(binraster) / (dt * 1000);
mltp = binraster;
mltp(~isnan(mltp)) = 1;
figure
A = axes;
hold on
for k = 1:tno   % convolve trial-wise
    spks = find(binraster(k,:)==1);
    spno = length(spks);
    for t = 1:spno
        spi = spks(t);
        tbinraster = zeros(1,tl);
        tbinraster(spi) = 1;
        if prob(spi) > 1
            keyboard
        end
        wbh = gausswin(9,prob(spi)*50);
        wbh = wbh / sum(wbh);
        plot(wbh,'Color',rand(1,3))
        
        agvd(k,:) = agvd(k,:) + filtfilt(wbh,1,tbinraster);
    end
end
psth = nanmean(agvd.*mltp) / dt;
psth_err = nanstd(agvd.*mltp/dt) / sqrt(tno);

% Plot
H1 = figure;
subplot(211)
imagesc(agvd);
% set(gca,'XTick',[])
subplot(212)
errorshade(time,psth,psth_err,'LineColor','k','LineWidth',2,'ShadeColor','grey')
xlim([time(1) time(end)])   
axes('Position',[0.7 0.35 0.2 0.2])
Ls = findobj(allchild(A),'type','line');   % copyobj runs to Matlab bug
plot(cell2mat(get(Ls,'XData'))',cell2mat(get(Ls,'YData'))')