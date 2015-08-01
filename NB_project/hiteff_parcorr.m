function hiteff_parcorr
%HITEFF_PARCORR
%   HITEFF_PARCOR calculates partial correlations between recording depth
%   (dorso-ventral coordinate) and reward-elicited response (see
%   DEPTH_VS_HITRESPONSE) while controlling for training history (see
%   TRAINING_VS_HITRESPONSE).
%
%   See also PARCOR, DEPTH_VS_HITRESPONSE and TRAINING_VS_HITRESPONSE.

% Load variables
global DATAPATH
load([DATAPATH 'NB\depth_vs_hitresponse_newdata\parcorr\depthvars.mat'])   % depth
allChAT1 = allChAT;
load([DATAPATH 'NB\depth_vs_hitresponse_newdata\parcorr\trainingvars.mat'])   % training history
if ~isequal(allChAT1,allChAT)
    error('hiteff_parcorr:inputArgMismatch','Variables should contain the cells in corresponding order.')
end
Hiteff2 = RelEff;

% Partail correlation
partialcorr(depth,Hiteff2',daynum)
parcor(depth,Hiteff2,daynum)
parcor(depth,Hiteff2,sessionnum)
parcor(depth,Hiteff2,trialnum)
parcor(depth,Hiteff2,fanum)
parcor(depth,Hiteff2,hitnum)

parcor(daynum,Hiteff2,depth)