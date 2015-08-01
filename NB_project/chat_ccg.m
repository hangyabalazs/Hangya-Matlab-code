function chat_ccg
%CHAT_CCG   Auto- and cross-correlations.
%   CHAT_CCG calculated and saves auto- and cross-correlations for
%   choinergic and putative cholinergic cells. See NBACG and NBCCG for
%   details.
%
%   See also NBACG and NBCCG.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013

% Cholinergic neurons
ChAT1 = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
pChAT1 = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT1 pChAT1];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Auto-correlation
% acg(allChAT,'issave',true)
% keyboard

% Cholinergic neurons
ChAT2 = selectcell(['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
pChAT2 = selectcell(['"pChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [allChAT ChAT2 pChAT2];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Cross-correlation
ccg(allChAT,'issave',true,'whichcells','allpairs')
keyboard