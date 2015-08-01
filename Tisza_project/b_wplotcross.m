function b_wplotcross(wave_cross,f)
%WPLOTCROSS    Plots crosswavelet power with pointwise maximums.
%   WPLOTCROSS(W,F) plots power of W wavelet, calculates and visualizes
%   pointwise maximum locations of power and transformates y-axis from
%   scales to frequency using F scale vector.
%
%   See also WCROSSWAVELET and RESCALEAXIS.

% Power and phase
wave_power = abs(wave_cross) .^ 2;
wave_phase = angle(wave_cross);
sw2 = size(wave_power,2);

% Maximum localizations
maxes = max(wave_power);
maxloc = zeros(1,sw2);
for i = 1:sw2
    maxloc(i) = find(wave_power(:,i)==maxes(i));
end

% Plot
H = figure;
imagesc(wave_power)
hold on
plot(maxloc,'Color','m','Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 1 1],...
    'MarkerSize',2,'LineStyle','none')
b_rescaleaxis('Y',f)