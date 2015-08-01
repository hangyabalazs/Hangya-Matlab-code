function varargout = czspikesort(varargin)
% CZSPIKESORT M-file for czspikesort.fig
%      CZSPIKESORT, by itself, creates a new CZSPIKESORT or raises the existing
%      singleton*.
%
%      H = CZSPIKESORT returns the handle to a new CZSPIKESORT or the handle to
%      the existing singleton*.
%
%      CZSPIKESORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CZSPIKESORT.M with the given input arguments.
%
%      CZSPIKESORT('Property','Value',...) creates a new CZSPIKESORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before czspikesort_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to czspikesort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%   CZSPIKESORT provides a graphical user interface to manually separate
%   spike clusters from tetrode recordings. It plots 2D views of the
%   clusters, auto- and cross-correlations and spike waveforms. It
%   calculates Isolation Distance and L-ratio for all clusters.
%       'Create' menu - calculates and saves PCA coefficients
%       'Load' menu - loads PCA results
%       'Import' menu - imports previously exported clusters
%       'Export' menu - exports clusters
%       'Extended export' menu - exports clusters, auto- and
%           crosscorrelations, spatial correlations, transmission success 
%           rates and phase histograms
%       'Place field' menu - plots rate map
%   czspikesort_OpeningFcn should be adjusted to the imported file.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help czspikesort

% Last Modified by GUIDE v2.5 18-Feb-2010 14:09:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @czspikesort_OpeningFcn, ...
                   'gui_OutputFcn',  @czspikesort_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% -------------------------------------------------------------------------
function czspikesort_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for czspikesort
handles.output = hObject;

% File names
handles.Nttfile = 'X:\In_Vivo\balazs\_analysis\Czurko\Ntt\Hux047_R_Sc3.Ntt';
handles.inpdir = 'X:\In_Vivo\balazs\_analysis\Czurko\discriminated2\neg\hux047-day06-tr1-01_o\';
handles.filename = 'hux047-day06-tr1-01_o';
handles.corrfile = [handles.inpdir 'nr_02_1.mat'];

% Set popup menu string
str = [{'Red'}; {'Green'}; {'Blue'}; {'Magenta'}; {'Yellow'}];
set(handles.popupmenu_activecell,'String',str)

% Initialize feature names
feature_field_names = [{'EL1 PC1'} {'el1_pc1'}; {'EL2 PC1'} {'el2_pc1'}; ...
    {'EL3 PC1'} {'el3_pc1'}; {'EL4 PC1'} {'el4_pc1'}];
feature_field_names = [feature_field_names; {'EL1 PC2'} {'el1_pc2'}; {'EL2 PC2'} {'el2_pc2'}; ...
    {'EL3 PC2'} {'el3_pc2'}; {'EL4 PC2'} {'el4_pc2'}];
feature_field_names = [feature_field_names; {'EL1 PC3'} {'el1_pc3'}; {'EL2 PC3'} {'el2_pc3'}; ...
    {'EL3 PC3'} {'el3_pc3'}; {'EL4 PC3'} {'el4_pc3'}];
feature_field_names = [feature_field_names; {'EL1 PC4'} {'el1_pc4'}; {'EL2 PC4'} {'el2_pc4'}; ...
    {'EL3 PC4'} {'el3_pc4'}; {'EL4 PC4'} {'el4_pc4'}];
feature_field_names = [feature_field_names; {'EL1 PCE1'} {'el1_pce1'}; {'EL2 PCE1'} {'el2_pce1'}; ...
    {'EL3 PCE1'} {'el3_pce1'}; {'EL4 PCE1'} {'el4_pce1'}];
feature_field_names = [feature_field_names; {'EL1 PCE2'} {'el1_pce2'}; {'EL2 PCE2'} {'el2_pce2'}; ...
    {'EL3 PCE2'} {'el3_pce2'}; {'EL4 PCE2'} {'el4_pce2'}];
feature_field_names = [feature_field_names; {'EL1 PCE3'} {'el1_pce3'}; {'EL2 PCE3'} {'el2_pce3'}; ...
    {'EL3 PCE3'} {'el3_pce3'}; {'EL4 PCE3'} {'el4_pce3'}];
feature_field_names = [feature_field_names; {'EL1 PCE4'} {'el1_pce4'}; {'EL2 PCE4'} {'el2_pce4'}; ...
    {'EL3 PCE4'} {'el3_pce4'}; {'EL4 PCE4'} {'el4_pce4'}];
feature_field_names = [feature_field_names; {'EL1 PCA1'} {'el1_pca1'}; {'EL2 PCA1'} {'el2_pca1'}; ...
    {'EL3 PCA1'} {'el3_pca1'}; {'EL4 PCA1'} {'el4_pca1'}];
feature_field_names = [feature_field_names; {'EL1 PCA2'} {'el1_pca2'}; {'EL2 PCA2'} {'el2_pca2'}; ...
    {'EL3 PCA2'} {'el3_pca2'}; {'EL4 PCA2'} {'el4_pca2'}];
feature_field_names = [feature_field_names; {'EL1 PCA3'} {'el1_pca3'}; {'EL2 PCA3'} {'el2_pca3'}; ...
    {'EL3 PCA3'} {'el3_pca3'}; {'EL4 PCA3'} {'el4_pca3'}];
feature_field_names = [feature_field_names; {'EL1 PCA4'} {'el1_pca4'}; {'EL2 PCA4'} {'el2_pca4'}; ...
    {'EL3 PCA4'} {'el3_pca4'}; {'EL4 PCA4'} {'el4_pca4'}];
feature_field_names = [feature_field_names; {'EL1 PEAK'} {'el1_peak'}; {'EL2 PEAK'} {'el2_peak'}; ...
    {'EL3 PEAK'} {'el3_peak'}; {'EL4 PEAK'} {'el4_peak'}];
feature_field_names = [feature_field_names; {'EL1 VALLEY'} {'el1_valley'}; {'EL2 VALLEY'} {'el2_valley'}; ...
    {'EL3 VALLEY'} {'el3_valley'}; {'EL4 VALLEY'} {'el4_valley'}];
feature_field_names = [feature_field_names; {'EL1 AMP'} {'el1_amp'}; {'EL2 AMP'} {'el2_amp'}; ...
    {'EL3 AMP'} {'el3_amp'}; {'EL4 AMP'} {'el4_amp'}];
feature_field_names = [feature_field_names; {'EL1 E'} {'el1_energy'}; {'EL2 E'} {'el2_energy'}; ...
    {'EL3 E'} {'el3_energy'}; {'EL4 E'} {'el4_energy'}];
feature_field_names = [feature_field_names; {'EL1 A'} {'el1_area'}; {'EL2 A'} {'el2_area'}; ...
    {'EL3 A'} {'el3_area'}; {'EL4 A'} {'el4_area'}];
feature_field_names = [feature_field_names; {'EL1 P1'} {'el1_params1'}; {'EL2 P1'} {'el2_params1'}; ...
    {'EL3 P1'} {'el3_params1'}; {'EL4 P1'} {'el4_params1'}];
feature_field_names = [feature_field_names; {'EL1 P2'} {'el1_params2'}; {'EL2 P2'} {'el2_params2'}; ...
    {'EL3 P2'} {'el3_params2'}; {'EL4 P2'} {'el4_params2'}];
feature_field_names = [feature_field_names; {'EL1 SLICE'} {'el1_slice'}; {'EL2 SLICE'} {'el2_slice'}; ...
    {'EL3 SLICE'} {'el3_slice'}; {'EL4 SLICE'} {'el4_slice'}];
handles.feature_field_names = feature_field_names;
handles.xfeature = 'EL1 PC1';
handles.yfeatire = 'EL2 PC1';

% Initialize cell codes
handles.cell_codes = [1 2 9 13];
handles.pair_codes = [1 2; 1 9; 1 13];
handles.plot_codes = [{'r.'} {'g.'} {'b.'} {'m.'} {'y.'}];

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------
function varargout = czspikesort_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
% OPEN
% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)

% Add Nlx converter to Matlab path
dr = 'C:\MATLAB_R2007a\work\Nlx_converter';
addpath(genpath(dr))

% Import
[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints NlxHeader] = ...
    Nlx2MatSpike(handles.Nttfile,[1 1 1 1 1],1,1,1);

% Set variables
handles.TimeStamps = TimeStamps;
handles.CellNumbers = CellNumbers;
handles.Params = Params;
handles.DataPoints = DataPoints;
handles.CellNumbers0 = CellNumbers;
guidata(hObject,handles)

% PCA
DP1 = squeeze(DataPoints(:,1,:))';
[coeff_DP1, scores_DP1] = princomp(DP1);
DP2 = squeeze(DataPoints(:,2,:))';
[coeff_DP2, scores_DP2] = princomp(DP2);
DP3 = squeeze(DataPoints(:,3,:))';
[coeff_DP3, scores_DP3] = princomp(DP3);
DP4 = squeeze(DataPoints(:,4,:))';
[coeff_DP4, scores_DP4] = princomp(DP4);

% Slice
n = size(DataPoints,3);
Slice = zeros(4,n);
for z = 1:4
    slope = zeros(1,n);
    slope2 = zeros(1,n);
    intercept = zeros(1,n);
    for k = 1:n
        wn = squeeze(DataPoints(:,z,k))';
        mxinx = find(wn==max(wn));
        mxinx = mxinx(end);
        linlen = 5;
        lm = min(12,length(wn)-mxinx-linlen+1);
        lsrate = 10000;
        time = linspace(0,linlen/lsrate,linlen+1);
        r = zeros(1,lm);
        gr = zeros(1,lm);
        icp = zeros(1,lm);
        for m = 1:lm
            [gr0,icp0,err0] = linefit(time,wn(mxinx+m-1:mxinx+m+linlen-1));
            r(m) = err0;
            gr(m) = gr0;
            icp(m) = icp0;
        end
        r(gr>0) = 1000;
        lsr = min(r);
        grinx = find(r==lsr);
%         lsr = min(gr);
%         grinx = find(gr==lsr);
        if isempty(grinx);
            slope(k) = 500;
            slope2(k) = 0;
            intercept(k) = 0;
        else
            grinx = grinx(end);
            slope(k) = gr(grinx) / 1000;   % slope is given in mV/ms
            slope2(k) = gr(grinx);
            intercept(k) = icp(grinx);
        end

%         figure
%         plot(wn)
%         hold on
%         t = [mxinx+grinx-1:mxinx-1+grinx+linlen];
%         y = time .* slope2(k) + intercept(k);
%         plot(t,y,'r')
%         min(r)
    end
    Slice(z,:) = slope;
end

% Peak, valley
Param1 = squeeze(max(DataPoints,[],1));
Param2 = squeeze(min(DataPoints,[],1));
Param = [Param1; Param2];

% Energy
l1 = squeeze(DataPoints(:,1,:));
E1 = sum(l1.^2,1) / size(l1,2);
l2 = squeeze(DataPoints(:,2,:));
E2 = sum(l2.^2,1) / size(l2,2);
l3 = squeeze(DataPoints(:,3,:));
E3 = sum(l3.^2,1) / size(l3,2);
l4 = squeeze(DataPoints(:,4,:));
E4 = sum(l4.^2,1) / size(l4,2);

% Area
l1 = squeeze(DataPoints(:,1,:));
A1 = sum(abs(l1),1);
l2 = squeeze(DataPoints(:,2,:));
A2 = sum(abs(l2),1);
l3 = squeeze(DataPoints(:,3,:));
A3 = sum(abs(l3),1);
l4 = squeeze(DataPoints(:,4,:));
A4 = sum(abs(l4),1);

% PCA of energy normalized waveforms
DPE1 = squeeze(DataPoints(:,1,:))' ./ repmat(E1',1,32);
[coeff_DPE1, scores_DPE1] = princomp(DPE1);
DPE2 = squeeze(DataPoints(:,2,:))' ./ repmat(E2',1,32);
[coeff_DPE2, scores_DPE2] = princomp(DPE2);
DPE3 = squeeze(DataPoints(:,3,:))' ./ repmat(E3',1,32);
[coeff_DPE3, scores_DPE3] = princomp(DPE3);
DPE4 = squeeze(DataPoints(:,4,:))' ./ repmat(E4',1,32);
[coeff_DPE4, scores_DPE4] = princomp(DPE4);

% PCA of area normalized waveforms
DPA1 = squeeze(DataPoints(:,1,:))' ./ repmat(A1',1,32);
[coeff_DPA1, scores_DPA1] = princomp(DPA1);
DPA2 = squeeze(DataPoints(:,2,:))' ./ repmat(A2',1,32);
[coeff_DPA2, scores_DPA2] = princomp(DPA2);
DPA3 = squeeze(DataPoints(:,3,:))' ./ repmat(A3',1,32);
[coeff_DPA3, scores_DPA3] = princomp(DPA3);
DPA4 = squeeze(DataPoints(:,4,:))' ./ repmat(A4',1,32);
[coeff_DPA4, scores_DPA4] = princomp(DPA4);

% Plot
S1 = scores_DP1(:,1);
S2 = scores_DP3(:,1);
axes(handles.axes1)
ach = allchild(gca);
delete(ach)     % clear axes
plot(S1,S2,'k.','MarkerSize',4)
hold on

for k = 1:length(handles.cell_codes)
    inx = find(CellNumbers==handles.cell_codes(k));
    plot(S1(inx),S2(inx),handles.plot_codes{k},'MarkerSize',4)
end
handles.featS1 = S1;
handles.featS2 = S2;

% Pass variables
handles.x_feature = 0;
handles.y_feature = 0;
handles.el1_pc1 = scores_DP1(:,1);
handles.el1_pc2 = scores_DP1(:,2);
handles.el1_pc3 = scores_DP1(:,3);
handles.el1_pc4 = scores_DP1(:,4);
handles.el2_pc1 = scores_DP2(:,1);
handles.el2_pc2 = scores_DP2(:,2);
handles.el2_pc3 = scores_DP2(:,3);
handles.el2_pc4 = scores_DP2(:,4);
handles.el3_pc1 = scores_DP3(:,1);
handles.el3_pc2 = scores_DP3(:,2);
handles.el3_pc3 = scores_DP3(:,3);
handles.el3_pc4 = scores_DP3(:,4);
handles.el4_pc1 = scores_DP4(:,1);
handles.el4_pc2 = scores_DP4(:,2);
handles.el4_pc3 = scores_DP4(:,3);
handles.el4_pc4 = scores_DP4(:,4);

handles.el1_pce1 = scores_DPE1(:,1);
handles.el1_pce2 = scores_DPE1(:,2);
handles.el1_pce3 = scores_DPE1(:,3);
handles.el1_pce4 = scores_DPE1(:,4);
handles.el2_pce1 = scores_DPE2(:,1);
handles.el2_pce2 = scores_DPE2(:,2);
handles.el2_pce3 = scores_DPE2(:,3);
handles.el2_pce4 = scores_DPE2(:,4);
handles.el3_pce1 = scores_DPE3(:,1);
handles.el3_pce2 = scores_DPE3(:,2);
handles.el3_pce3 = scores_DPE3(:,3);
handles.el3_pce4 = scores_DPE3(:,4);
handles.el4_pce1 = scores_DPE4(:,1);
handles.el4_pce2 = scores_DPE4(:,2);
handles.el4_pce3 = scores_DPE4(:,3);
handles.el4_pce4 = scores_DPE4(:,4);

handles.el1_pca1 = scores_DPA1(:,1);
handles.el1_pca2 = scores_DPA1(:,2);
handles.el1_pca3 = scores_DPA1(:,3);
handles.el1_pca4 = scores_DPA1(:,4);
handles.el2_pca1 = scores_DPA2(:,1);
handles.el2_pca2 = scores_DPA2(:,2);
handles.el2_pca3 = scores_DPA2(:,3);
handles.el2_pca4 = scores_DPA2(:,4);
handles.el3_pca1 = scores_DPA3(:,1);
handles.el3_pca2 = scores_DPA3(:,2);
handles.el3_pca3 = scores_DPA3(:,3);
handles.el3_pca4 = scores_DPA3(:,4);
handles.el4_pca1 = scores_DPA4(:,1);
handles.el4_pca2 = scores_DPA4(:,2);
handles.el4_pca3 = scores_DPA4(:,3);
handles.el4_pca4 = scores_DPA4(:,4);

handles.el1_peak = Param(1,:);
handles.el2_peak = Param(2,:);
handles.el3_peak = Param(3,:);
handles.el4_peak = Param(4,:);
handles.el1_valley = Param(5,:);
handles.el2_valley = Param(6,:);
handles.el3_valley = Param(7,:);
handles.el4_valley = Param(8,:);
handles.el1_amp = Param(1,:) - Param(5,:);
handles.el2_amp = Param(2,:) - Param(6,:);
handles.el3_amp = Param(3,:) - Param(7,:);
handles.el4_amp = Param(4,:) - Param(8,:);

handles.el1_energy = E1;
handles.el2_energy = E2;
handles.el3_energy = E3;
handles.el4_energy = E4;
handles.el1_area = A1;
handles.el2_area = A2;
handles.el3_area = A3;
handles.el4_area = A4;
handles.el1_params1 = Params(1,:);
handles.el2_params1 = Params(2,:);
handles.el3_params1 = Params(3,:);
handles.el4_params1 = Params(4,:);
handles.el1_params2 = Params(5,:);
handles.el2_params2 = Params(6,:);
handles.el3_params2 = Params(7,:);
handles.el4_params2 = Params(8,:);

handles.el1_slice = Slice(1,:);
handles.el2_slice = Slice(2,:);
handles.el3_slice = Slice(3,:);
handles.el4_slice = Slice(4,:);
guidata(handles.figure1,handles)

% Enable feature space setting
set(handles.sorting,'Enable','on')
set(handles.analysis,'Enable','on')
set(handles.export,'Enable','on')
set(handles.extended_export,'Enable','on')

% Assign global GUI_HANDLE
global GUI_HANDLE
GUI_HANDLE = handles.figure1;

% Set default features
feature_names = {};
feature_vector = struct();
feature_vector.pc1 = 1;
feature_names = [feature_names {'EL1 PC1'} {'EL2 PC1'} {'EL3 PC1'} {'EL4 PC1'}];
feature_vector.pc2 = 1;
feature_names = [feature_names {'EL1 PC2'} {'EL2 PC2'} {'EL3 PC2'} {'EL4 PC2'}];
feature_vector.pc3 = 0;
feature_vector.pc4 = 0;
feature_vector.pce1 = 0;
feature_vector.pce2 = 0;
feature_vector.pce3 = 0;
feature_vector.pce4 = 0;
feature_vector.pca1 = 0;
feature_vector.pca2 = 0;
feature_vector.pca3 = 0;
feature_vector.pca4 = 0;
feature_vector.peak = 1;
feature_names = [feature_names {'EL1 PEAK'} {'EL2 PEAK'} {'EL3 PEAK'} {'EL4 PEAK'}];
feature_vector.valley = 1;
feature_names = [feature_names {'EL1 VALLEY'} {'EL2 VALLEY'} {'EL3 VALLEY'} {'EL4 VALLEY'}];
feature_vector.amp = 0;
feature_vector.energy = 0;
feature_vector.area = 0;
feature_vector.params = 0;
global FEATURE_NAMES
FEATURE_NAMES = feature_names;
global FEATURE_VECTOR
FEATURE_VECTOR = feature_vector;

% Set buttondown function
bdf = 'czspikesort(''Axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)

% Set keypress function
kpf = 'czspikesort(''nextfeature'')';
set(handles.figure1,'KeyPressFcn',kpf)
guidata(handles.figure1,handles)

% Set active cell
popupmenu_activecell_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function create_Callback(hObject, eventdata, handles)

% Add Nlx converter to Matlab path
dr = 'C:\MATLAB_R2007a\work\Nlx_converter';
addpath(genpath(dr))

% Import
[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints NlxHeader] = ...
    Nlx2MatSpike(handles.Nttfile,[1 1 1 1 1],1,1,1);
inpdir = handles.inpdir;

% Timestamp correction
load(handles.corrfile)
data = data';
inx = find(CellNumbers==1);
TimeStamps2 = TimeStamps - TimeStamps(1);
data2 = TimeStamps2(inx);
min(data-data2/1000000)
max(data-data2/1000000)
corrv = mean(data-data2/1000000)
handles.corrv = corrv;

% PCA
DP1 = squeeze(DataPoints(:,1,:))';
[coeff_DP1, scores_DP1] = princomp(DP1);
DP2 = squeeze(DataPoints(:,2,:))';
[coeff_DP2, scores_DP2] = princomp(DP2);
DP3 = squeeze(DataPoints(:,3,:))';
[coeff_DP3, scores_DP3] = princomp(DP3);
DP4 = squeeze(DataPoints(:,4,:))';
[coeff_DP4, scores_DP4] = princomp(DP4);

% Slice
n = size(DataPoints,3);
Slice = zeros(4,n);
for z = 1:4
    slope = zeros(1,n);
    slope2 = zeros(1,n);
    intercept = zeros(1,n);
    for k = 1:n
        wn = squeeze(DataPoints(:,z,k))';
        mxinx = find(wn==max(wn));
        mxinx = mxinx(end);
        linlen = 5;
        lm = min(12,length(wn)-mxinx-linlen+1);
        lsrate = 10000;
        time = linspace(0,linlen/lsrate,linlen+1);
        r = zeros(1,lm);
        gr = zeros(1,lm);
        icp = zeros(1,lm);
        for m = 1:lm
            [gr0,icp0,err0] = linefit(time,wn(mxinx+m-1:mxinx+m+linlen-1));
            r(m) = err0;
            gr(m) = gr0;
            icp(m) = icp0;
        end
        r(gr>0) = 1000;
        lsr = min(r);
        grinx = find(r==lsr);
%         lsr = min(gr);
%         grinx = find(gr==lsr);
        if isempty(grinx);
            slope(k) = 500;
            slope2(k) = 0;
            intercept(k) = 0;
        else
            grinx = grinx(end);
            slope(k) = gr(grinx) / 1000;   % slope is given in mV/ms
            slope2(k) = gr(grinx);
            intercept(k) = icp(grinx);
        end

%         figure
%         plot(wn)
%         hold on
%         t = [mxinx+grinx-1:mxinx-1+grinx+linlen];
%         y = time .* slope2(k) + intercept(k);
%         plot(t,y,'r')
%         min(r)
    end
    Slice(z,:) = slope;
end
Slice(Slice<-10000) = -10000;
Slice(Slice>500) = 500;

% Peak, valley
Param1 = squeeze(max(DataPoints,[],1));
Param2 = squeeze(min(DataPoints,[],1));
Param = [Param1; Param2];

% Energy
l1 = squeeze(DataPoints(:,1,:));
E1 = sum(l1.^2,1) / size(l1,2);
l2 = squeeze(DataPoints(:,2,:));
E2 = sum(l2.^2,1) / size(l2,2);
l3 = squeeze(DataPoints(:,3,:));
E3 = sum(l3.^2,1) / size(l3,2);
l4 = squeeze(DataPoints(:,4,:));
E4 = sum(l4.^2,1) / size(l4,2);

% Area
l1 = squeeze(DataPoints(:,1,:));
A1 = sum(abs(l1),1);
l2 = squeeze(DataPoints(:,2,:));
A2 = sum(abs(l2),1);
l3 = squeeze(DataPoints(:,3,:));
A3 = sum(abs(l3),1);
l4 = squeeze(DataPoints(:,4,:));
A4 = sum(abs(l4),1);

% PCA of energy normalized waveforms
DPE1 = squeeze(DataPoints(:,1,:))' ./ repmat(E1',1,32);
[coeff_DPE1, scores_DPE1] = princomp(DPE1);
DPE2 = squeeze(DataPoints(:,2,:))' ./ repmat(E2',1,32);
[coeff_DPE2, scores_DPE2] = princomp(DPE2);
DPE3 = squeeze(DataPoints(:,3,:))' ./ repmat(E3',1,32);
[coeff_DPE3, scores_DPE3] = princomp(DPE3);
DPE4 = squeeze(DataPoints(:,4,:))' ./ repmat(E4',1,32);
[coeff_DPE4, scores_DPE4] = princomp(DPE4);

% PCA of area normalized waveforms
DPA1 = squeeze(DataPoints(:,1,:))' ./ repmat(A1',1,32);
[coeff_DPA1, scores_DPA1] = princomp(DPA1);
DPA2 = squeeze(DataPoints(:,2,:))' ./ repmat(A2',1,32);
[coeff_DPA2, scores_DPA2] = princomp(DPA2);
DPA3 = squeeze(DataPoints(:,3,:))' ./ repmat(A3',1,32);
[coeff_DPA3, scores_DPA3] = princomp(DPA3);
DPA4 = squeeze(DataPoints(:,4,:))' ./ repmat(A4',1,32);
[coeff_DPA4, scores_DPA4] = princomp(DPA4);

% Save
resdir = [inpdir 'Recluster0\'];     % result directory
if ~isdir(resdir)
    mkdir(resdir)
end
fn = [resdir 'features.mat'];       % result file
save(fn,'scores_DP1','scores_DP2','scores_DP3','scores_DP4',...
    'scores_DPA1','scores_DPA2','scores_DPA3','scores_DPA4',...
    'scores_DPE1','scores_DPE2','scores_DPE3','scores_DPE4',...
    'A1','A2','A3','A4','E1','E2','E3','E4','Param','Params','Slice')
disp('Ready.')

% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)

% Load
inpdir = handles.inpdir;
resdir = [inpdir 'Recluster0_sc3\'];     % load directory
fn = [resdir 'features.mat'];       % result file
load(fn)

% Add Nlx converter to Matlab path
dr = 'C:\MATLAB_R2007a\work\Nlx_converter';
addpath(genpath(dr))

% Import
[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints NlxHeader] = ...
    Nlx2MatSpike(handles.Nttfile,...
    [1 1 1 1 1],1,1,1);


% Set variables
handles.TimeStamps = TimeStamps;
handles.CellNumbers = CellNumbers;
handles.Params = Params;
handles.DataPoints = DataPoints;
handles.CellNumbers0 = CellNumbers;
guidata(hObject,handles)

% Timestamp correction
load(handles.corrfile)
data = data';
inx = find(handles.CellNumbers==1);
TimeStamps2 = TimeStamps - TimeStamps(1);
data2 = TimeStamps2(inx);
min(data-data2/1000000)
max(data-data2/1000000)
handles.corrv = mean(data-data2/1000000);

% Plot
S1 = scores_DP1(:,1);
S2 = scores_DP3(:,1);
axes(handles.axes1)
ach = allchild(gca);
delete(ach)     % clear axes
plot(S1,S2,'k.','MarkerSize',4)
hold on

for k = 1:length(handles.cell_codes)
    inx = find(CellNumbers==handles.cell_codes(k));
    plot(S1(inx),S2(inx),handles.plot_codes{k},'MarkerSize',4)
end
handles.featS1 = S1;
handles.featS2 = S2;

% Pass variables
handles.x_feature = 0;
handles.y_feature = 0;
handles.el1_pc1 = scores_DP1(:,1);
handles.el1_pc2 = scores_DP1(:,2);
handles.el1_pc3 = scores_DP1(:,3);
handles.el1_pc4 = scores_DP1(:,4);
handles.el2_pc1 = scores_DP2(:,1);
handles.el2_pc2 = scores_DP2(:,2);
handles.el2_pc3 = scores_DP2(:,3);
handles.el2_pc4 = scores_DP2(:,4);
handles.el3_pc1 = scores_DP3(:,1);
handles.el3_pc2 = scores_DP3(:,2);
handles.el3_pc3 = scores_DP3(:,3);
handles.el3_pc4 = scores_DP3(:,4);
handles.el4_pc1 = scores_DP4(:,1);
handles.el4_pc2 = scores_DP4(:,2);
handles.el4_pc3 = scores_DP4(:,3);
handles.el4_pc4 = scores_DP4(:,4);

handles.el1_pce1 = scores_DPE1(:,1);
handles.el1_pce2 = scores_DPE1(:,2);
handles.el1_pce3 = scores_DPE1(:,3);
handles.el1_pce4 = scores_DPE1(:,4);
handles.el2_pce1 = scores_DPE2(:,1);
handles.el2_pce2 = scores_DPE2(:,2);
handles.el2_pce3 = scores_DPE2(:,3);
handles.el2_pce4 = scores_DPE2(:,4);
handles.el3_pce1 = scores_DPE3(:,1);
handles.el3_pce2 = scores_DPE3(:,2);
handles.el3_pce3 = scores_DPE3(:,3);
handles.el3_pce4 = scores_DPE3(:,4);
handles.el4_pce1 = scores_DPE4(:,1);
handles.el4_pce2 = scores_DPE4(:,2);
handles.el4_pce3 = scores_DPE4(:,3);
handles.el4_pce4 = scores_DPE4(:,4);

handles.el1_pca1 = scores_DPA1(:,1);
handles.el1_pca2 = scores_DPA1(:,2);
handles.el1_pca3 = scores_DPA1(:,3);
handles.el1_pca4 = scores_DPA1(:,4);
handles.el2_pca1 = scores_DPA2(:,1);
handles.el2_pca2 = scores_DPA2(:,2);
handles.el2_pca3 = scores_DPA2(:,3);
handles.el2_pca4 = scores_DPA2(:,4);
handles.el3_pca1 = scores_DPA3(:,1);
handles.el3_pca2 = scores_DPA3(:,2);
handles.el3_pca3 = scores_DPA3(:,3);
handles.el3_pca4 = scores_DPA3(:,4);
handles.el4_pca1 = scores_DPA4(:,1);
handles.el4_pca2 = scores_DPA4(:,2);
handles.el4_pca3 = scores_DPA4(:,3);
handles.el4_pca4 = scores_DPA4(:,4);

handles.el1_peak = Param(1,:);
handles.el2_peak = Param(2,:);
handles.el3_peak = Param(3,:);
handles.el4_peak = Param(4,:);
handles.el1_valley = Param(5,:);
handles.el2_valley = Param(6,:);
handles.el3_valley = Param(7,:);
handles.el4_valley = Param(8,:);
handles.el1_amp = Param(1,:) - Param(5,:);
handles.el2_amp = Param(2,:) - Param(6,:);
handles.el3_amp = Param(3,:) - Param(7,:);
handles.el4_amp = Param(4,:) - Param(8,:);

handles.el1_energy = E1;
handles.el2_energy = E2;
handles.el3_energy = E3;
handles.el4_energy = E4;
handles.el1_area = A1;
handles.el2_area = A2;
handles.el3_area = A3;
handles.el4_area = A4;
handles.el1_params1 = Params(1,:);
handles.el2_params1 = Params(2,:);
handles.el3_params1 = Params(3,:);
handles.el4_params1 = Params(4,:);
handles.el1_params2 = Params(5,:);
handles.el2_params2 = Params(6,:);
handles.el3_params2 = Params(7,:);
handles.el4_params2 = Params(8,:);

handles.el1_slice = Slice(1,:);
handles.el2_slice = Slice(2,:);
handles.el3_slice = Slice(3,:);
handles.el4_slice = Slice(4,:);
guidata(handles.figure1,handles)

% Enable feature space setting
set(handles.sorting,'Enable','on')
set(handles.analysis,'Enable','on')
set(handles.export,'Enable','on')
set(handles.extended_export,'Enable','on')

% Assign global GUI_HANDLE
global GUI_HANDLE
GUI_HANDLE = handles.figure1;

% Set default features
feature_names = {};
feature_vector = struct();
feature_vector.pc1 = 1;
feature_names = [feature_names {'EL1 PC1'} {'EL2 PC1'} {'EL3 PC1'} {'EL4 PC1'}];
feature_vector.pc2 = 1;
feature_names = [feature_names {'EL1 PC2'} {'EL2 PC2'} {'EL3 PC2'} {'EL4 PC2'}];
feature_vector.pc3 = 0;
feature_vector.pc4 = 0;
feature_vector.pce1 = 0;
feature_vector.pce2 = 0;
feature_vector.pce3 = 0;
feature_vector.pce4 = 0;
feature_vector.pca1 = 0;
feature_vector.pca2 = 0;
feature_vector.pca3 = 0;
feature_vector.pca4 = 0;
feature_vector.peak = 1;
feature_names = [feature_names {'EL1 PEAK'} {'EL2 PEAK'} {'EL3 PEAK'} {'EL4 PEAK'}];
feature_vector.valley = 1;
feature_names = [feature_names {'EL1 VALLEY'} {'EL2 VALLEY'} {'EL3 VALLEY'} {'EL4 VALLEY'}];
feature_vector.amp = 0;
feature_vector.slice = 0;
feature_vector.energy = 0;
feature_vector.area = 0;
feature_vector.params = 0;
global FEATURE_NAMES
FEATURE_NAMES = feature_names;
global FEATURE_VECTOR
FEATURE_VECTOR = feature_vector;

% Set buttondown function
bdf = 'czspikesort(''Axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)

% Set keypress function
kpf = 'czspikesort(''nextfeature'')';
set(handles.figure1,'KeyPressFcn',kpf)
guidata(handles.figure1,handles)

% Set active cell
popupmenu_activecell_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function import_Callback(hObject, eventdata, handles)

% Add Nlx converter to Matlab path
dr = 'C:\MATLAB_R2007a\work\Nlx_converter';
addpath(genpath(dr))

% Import
impdir = [handles.inpdir 'Recluster_sc3\'];
fn = [impdir 'clusters.mat'];
load(fn)

% Assign variables
handles.DataPoints = DataPoints;
handles.CellNumbers = CellNumbers;
handles.TimeStamps = TimeStamps;
global FEATURE_NAMES
FEATURE_NAMES = FeatureNames;
global FEATURE_VECTOR
FEATURE_VECTOR = FeatureVector;
guidata(handles.figure1,handles)

% Clear axes
ach = allchild(handles.axes1);
delete(ach)

% Replot
replotsub(handles)

% Set active cell
popupmenu_activecell_Callback(hObject, eventdata, handles)

% Set buttondown function
bdf = 'czspikesort(''Axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function sorting_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function feature_space_Callback(hObject, eventdata, handles)

G = czfeaturespace;
uiwait(G)

% Get index of selected cell
cellinx = get(handles.popupmenu_activecell,'Value');

% Isolation Distance and L-ratio
[DM LrC] = isdist(cellinx,handles);
set(handles.edit_isolationdist,'String',num2str(DM))
set(handles.edit_Lratio,'String',num2str(LrC))

% --------------------------------------------------------------------
% KEYPRESS FUNCTION
% --------------------------------------------------------------------
function nextfeature

% Calling the appropriate functions
inp = get(gcf,'CurrentCharacter');
switch inp
    case 'y'    % y - switch y feature forward
        replot('y')
    case 'x'    % x - switch x feature forward
        replot('x')
    case 's'    % s - switch y feature backward
        replot('s')
    case 'd'    % d - switch x feature backward
        replot('d')
end

% --------------------------------------------------------------------
function replot(s)

% Get handles structure
global GUI_HANDLE
h = GUI_HANDLE;
handles = guidata(h);

% Load feature list
global FEATURE_NAMES
global FEATURE_VECTOR

% Reset feature axes
axes(handles.axes1)
ach = allchild(gca);
delete(ach)     % clear axes
switch s
    case 'x'
        handles.x_feature = handles.x_feature + 1;
    case 'y'
        handles.y_feature = handles.y_feature + 1;
    case 'd'
        handles.x_feature = handles.x_feature - 1;
    case 's'
        handles.y_feature = handles.y_feature - 1;
end
% if handles.x_feature == 0
%     handles.x_feature = 1;
% end
% if handles.y_feature == 0
%     handles.y_feature = 1;
% end
feature_num = length(FEATURE_NAMES);
handles.x_feature = mod2(handles.x_feature,feature_num);
handles.y_feature = mod2(handles.y_feature,feature_num);
handles.xfeature = FEATURE_NAMES{handles.x_feature};
handles.yfeature = FEATURE_NAMES{handles.y_feature};
guidata(h,handles)

% Replot
replotsub(handles)

% -------------------------------------------------------------------------
function replotsub(handles)

% Get feature variables
inx = find(strcmp(handles.xfeature,handles.feature_field_names(:,1)));
xfield = handles.feature_field_names{inx,2};
eval(['S1 = handles.' xfield ';']);
inx = find(strcmp(handles.yfeature,handles.feature_field_names(:,1)));
yfield = handles.feature_field_names{inx,2};
eval(['S2 = handles.' yfield ';']);

% Plot
axes(handles.axes1)
plot(S1,S2,'k.','MarkerSize',4)
hold on

for k = 1:length(handles.cell_codes)
    inx = find(handles.CellNumbers==handles.cell_codes(k));
    plot(S1(inx),S2(inx),handles.plot_codes{k},'MarkerSize',4)
end
xlabel(handles.xfeature)
ylabel(handles.yfeature)
handles.featS1 = S1;
handles.featS2 = S2;
guidata(handles.figure1,handles)

% -------------------------------------------------------------------------
% Popup menu for ACTIVE CELL SELECTION
% -------------------------------------------------------------------------
function popupmenu_activecell_Callback(hObject, eventdata, handles)

% Get index of selected cell
cellinx = get(handles.popupmenu_activecell,'Value');

% Plot waveforms
plotwaveforms(cellinx,handles)

% Plot auto-correlogram
plotauto(cellinx,handles)

% Cluster indeces
inx = find(handles.CellNumbers==handles.cell_codes(cellinx));

% Spike number
set(handles.edit_spikenumber,'String',num2str(length(inx)))

% Isolation Distance and L-ratio
[DM LrC] = isdist(cellinx,handles);
set(handles.edit_isolationdist,'String',num2str(DM))
set(handles.edit_Lratio,'String',num2str(LrC))

% -------------------------------------------------------------------------
function popupmenu_activecell_CreateFcn(hObject, eventdata, handles)

% Popupmenu controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function plotwaveforms(cellindex,handles)
% Plot waveforms.

% Cluster indeces
inx = find(handles.CellNumbers==handles.cell_codes(cellindex));
cst = ceil(length(inx)/500);

% Plot
DataPoints = handles.DataPoints;
axes(handles.axes2)
ach = allchild(gca);
delete(ach)     % clear axes
plot(squeeze(DataPoints(:,1,inx(1:cst:end))),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,1,inx),3)),'r')
plot(squeeze(quantile(DataPoints(:,1,inx),0.95,3)),'b')
plot(squeeze(quantile(DataPoints(:,1,inx),0.05,3)),'b')
xlim([0 32])

axes(handles.axes3)
ach = allchild(gca);
delete(ach)     % clear axes
plot(squeeze(DataPoints(:,2,inx(1:cst:end))),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,2,inx),3)),'r')
plot(squeeze(quantile(DataPoints(:,2,inx),0.95,3)),'b')
plot(squeeze(quantile(DataPoints(:,2,inx),0.05,3)),'b')
xlim([0 32])

axes(handles.axes4)
ach = allchild(gca);
delete(ach)     % clear axes
plot(squeeze(DataPoints(:,3,inx(1:cst:end))),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,3,inx),3)),'r')
plot(squeeze(quantile(DataPoints(:,3,inx),0.95,3)),'b')
plot(squeeze(quantile(DataPoints(:,3,inx),0.05,3)),'b')
xlim([0 32])

axes(handles.axes5)
ach = allchild(gca);
delete(ach)     % clear axes
plot(squeeze(DataPoints(:,4,inx(1:cst:end))),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,4,inx),3)),'r')
plot(squeeze(quantile(DataPoints(:,4,inx),0.95,3)),'b')
plot(squeeze(quantile(DataPoints(:,4,inx),0.05,3)),'b')
xlim([0 32])

% -------------------------------------------------------------------------
function plotauto(cellindex,handles)
% Plot auto-correlogram.

% Cluster indeces
inx = find(handles.CellNumbers==handles.cell_codes(cellindex));
vdisc = (handles.TimeStamps(inx) - handles.TimeStamps(1)) / 1000000;

% Plot
axes(handles.axes6)
ach = allchild(gca);
delete(ach)     % clear axes
lczacorr(vdisc)

% -------------------------------------------------------------------------
function lczacorr(ncc)
%CZACORR   Autocorrelation.
%   CZACORR(VD) calculates autocorrelogram for discriminated unit VD, using
%   a +-50 ms time window.
%
%   See also XCORR and CZXCORR.

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;

% Autocorrelogram
zunit1 = zeros(1,length(round(nc))+5);
rnc = round(nc);
rnc = rnc(rnc>0);
zunit1(rnc) = 1;
acr = xcorr(zunit1,0.05*sr);
acr(length(acr)/2+0.5) = [];

% Plot result
bar(linspace(-50,50,length(acr)),acr)
xlim([-50 50])

% -------------------------------------------------------------------------
function [H1 H2 trsc] = lczxcorr(ncc1,ncc2)
%CZXCORR   Crosscorrelation.
%   CZXCORR(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-50 ms time window. Crosscorrelogram is
%   normalized with a shuffled ISI crosscorrelogram to remove random
%   correlations. A significance level of p=0.0013 is indicated on the
%   result. Original and normalized crosscorrelograms are plotted.
%
%   [H1 H2] = CZXCORR(VD1,VD2) returns the handles of the resulting plots.
%
%   [H1 H2 TRSC] = CZXCORR(VD1,VD2) returns transmission success rate as
%   well. Note, that interneuron has to be the second input argument and 
%   pyramidal cell has to be the first input argument for transmission
%   success rate calculation.
%
%   See also XCORR and CZACORR.

% Input argument check
error(nargchk(2,2,nargin))

% Calculate spike times in milliseconds
sr = 1000;
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
nc1(nc1<0.5) = [];
nc2(nc2<0.5) = [];

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
ccr = xcorr(zunit2,zunit1,0.05*sr);     % 1->2; window: -50 ms - 50 ms
ccr = ccr / length(nc1);     % norm. with the no. of ref. events to get transmission prob.
H1 = figure;
bar(linspace(-50,50,length(ccr)),ccr)
set(gca,'XLim',[-50 50])

% ISI shuffle
% isi = diff(nc2);
% lit = length(isi);
% rp = randperm(lit);
% psi1 = [];
% for it = 1:lit
%     psi1 = [psi1 isi(rp(it))];
% end
% psvd = [nc2(1) nc2(1)+cumsum(psi1)];
% pzunit = zeros(1,length(round(psvd))+5);
% pzunit(round(psvd)) = 1;
for k = 1:10
    rnd = rand(1) * 100 + 50;
    shf = round(rnd);
    pzunit = [zunit2(shf+1:end) zunit2(1:shf)];
    
% Random crosscorrelogram
    pccr = xcorr(pzunit,zunit1,0.05*sr);
    pccr = pccr / length(nc1);
    str = ['p' num2str(k) '=pccr;'];
    eval(str)
end
pccr = mean([p1;p2;p3;p4;p5;p6;p7;p8;p9;p10]);
% figure;
% bar(pccr)

% Normalized crosscorrelogram
nccr = ccr - pccr;
trsc = max(nccr(45:55));      % transmission success rate
H2 = figure;
bar(linspace(-50,50,length(nccr)),nccr)
set(gca,'XLim',[-50 50])

mnc = mean(nccr);
sdc = std(nccr);
thr = mnc + 3 * sdc;
nthr = mnc - 3 * sdc;
line([-50 50],[thr thr],'Color','red')      % significance level: p=0.0013
line([-50 50],[nthr nthr],'Color','red')
line([-50 50],[nthr nthr],'Color','red')

% -------------------------------------------------------------------------
function [H1 H2] = lczxcorr2(ncc1,ncc2)
%CZXCORR2   Crosscorrelation.
%   CZXCORR(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-5 ms time window. Crosscorrelogram is
%   normalized with a shuffled ISI crosscorrelogram to remove random
%   correlations. A significance level of p=0.0013 is indicated on the
%   result. Original and normalized crosscorrelograms are plotted.
%
%   [H1 H2] = CZXCORR2(VD1,VD2) returns the handles of the resulting plots.
%
%   See also XCORR and CZACORR.

% Input argument check
error(nargchk(2,2,nargin))

% Calculate spike times in milliseconds
sr = 2000;  % because it is fed to ccr, it is fine to choose an sr other than 1000
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
nc1(nc1<0.5) = [];
nc2(nc2<0.5) = [];

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
ccr = xcorr(zunit2,zunit1,0.005*sr);     % 1->2; window: -5 ms - 5 ms
ccr = ccr / length(nc1);     % norm. with the no. of ref. events to get transmission prob.
H1 = figure;
bar(linspace(-5,5,length(ccr)),ccr)
set(gca,'XLim',[-5 5])

% ISI shuffle
% isi = diff(nc2);
% lit = length(isi);
% rp = randperm(lit);
% psi1 = [];
% for it = 1:lit
%     psi1 = [psi1 isi(rp(it))];
% end
% psvd = [nc2(1) nc2(1)+cumsum(psi1)];
% pzunit = zeros(1,length(round(psvd))+5);
% pzunit(round(psvd)) = 1;
for k = 1:10
    rnd = rand(1) * 100 + 50;
    shf = round(rnd);
    pzunit = [zunit2(shf+1:end) zunit2(1:shf)];
    
% Random crosscorrelogram
    pccr = xcorr(pzunit,zunit1,0.005*sr);
    pccr = pccr / length(nc1);
    str = ['p' num2str(k) '=pccr;'];
    eval(str)
end
pccr = mean([p1;p2;p3;p4;p5;p6;p7;p8;p9;p10]);
% figure;
% bar(pccr)

% Normalized crosscorrelogram
nccr = ccr - pccr;
H2 = figure;
bar(linspace(-5,5,length(nccr)),nccr)
set(gca,'XLim',[-5 5])

mnc = mean(nccr);
sdc = std(nccr);
thr = mnc + 3 * sdc;
nthr = mnc - 3 * sdc;
line([-5 5],[thr thr],'Color','red')      % significance level: p=0.0013
line([-5 5],[nthr nthr],'Color','red')
line([-5 5],[nthr nthr],'Color','red')

% -------------------------------------------------------------------------
function [DM LrC] = isdist(cellindex,handles)
% Isolation Distance.

% Cluster indeces
n = length(handles.CellNumbers);
inx = find(handles.CellNumbers==handles.cell_codes(cellindex));
cinx = find(handles.CellNumbers~=handles.cell_codes(cellindex));

% Feature matrix
global FEATURE_NAMES
X = [];
for k = 1:length(FEATURE_NAMES)
    featinx = find(strcmp(FEATURE_NAMES{k},handles.feature_field_names(:,1)));
    currfeat = handles.feature_field_names{featinx,2};
    eval(['currfeat_vector = handles.' currfeat ';']);
    if size(currfeat_vector,1) ~= 1
        currfeat_vector = currfeat_vector';
    end
    X = [X; currfeat_vector];
end

% Mean and covariance
XC = X(:,inx);
muC = mean(XC,2);
SC = cov(XC');
V = X - repmat(muC,1,n);
SCinv = inv(SC);

% Isolation Distance
D = mahal(X',XC');
DC = D(cinx);
sDC = sort(DC);
linx = length(inx);
if linx <= length(sDC)
    DM = sDC(length(inx));
else
    DM = 0;
end

% L-ratio
df = size(X,1);
LC = sum(1-chi2cdf(DC,df));
nc = size(XC,2);
LrC = LC / nc;

% -------------------------------------------------------------------------
function edit_spikenumber_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit_spikenumber_CreateFcn(hObject, eventdata, handles)

% Edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function edit_isolationdist_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit_isolationdist_CreateFcn(hObject, eventdata, handles)

% Edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function edit_Lratio_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit_Lratio_CreateFcn(hObject, eventdata, handles)

% Edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function radiobutton_replace_Callback(hObject, eventdata, handles)

set(handles.radiobutton_replace,'Value',1)
set(handles.radiobutton_add,'Value',0)
set(handles.radiobutton_subtract,'Value',0)

% -------------------------------------------------------------------------
function radiobutton_add_Callback(hObject, eventdata, handles)

set(handles.radiobutton_replace,'Value',0)
set(handles.radiobutton_add,'Value',1)
set(handles.radiobutton_subtract,'Value',0)

% -------------------------------------------------------------------------
function radiobutton_subtract_Callback(hObject, eventdata, handles)

set(handles.radiobutton_replace,'Value',0)
set(handles.radiobutton_add,'Value',0)
set(handles.radiobutton_subtract,'Value',1)


% -------------------------------------------------------------------------
% Buttondown function: RECLUSTER
% -------------------------------------------------------------------------
function Axes1_ButtonDownFcn(hObject, eventdata, handles)

% Delete lines
ach = allchild(handles.axes1);
lns = findobj(ach,'Color','c');
delete(lns)

% Set buttondown function
bdf = [];   % prevent from reexecuting while drwing the polygon
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)

% Draw polygon
seltyp = get(handles.figure1,'SelectionType');  % click type
if isequal(seltyp,'normal')
    point1 = get(handles.axes1,'CurrentPoint'); % button down detected
    point1x = point1(1,1);
    point1y = point1(1,2);
    point0x = point1(1,1);
    point0y = point1(1,2);
    L = [];
    bp = 1;
    while bp
        bp = waitforbuttonpress;    % wait for mouse click
    end
    seltyp2 = get(handles.figure1,'SelectionType');  % click type
    while isequal(seltyp2,'normal')
        point2 = get(handles.axes1,'CurrentPoint'); % button down detected
        point2x = point2(1,1);
        point2y = point2(1,2);
        L(end+1) = line([point1x point2x],[point1y point2y],'Color','c','LineWidth',2);
        point0x = [point0x point2x];
        point0y = [point0y point2y];
        point1x = point2x;
        point1y = point2y;
        bp = 1;
        while bp
            bp = waitforbuttonpress;    % wait for mouse click
        end
        seltyp2 = get(handles.figure1,'SelectionType');  % click type
    end
    L(end+1) = line([point2x point0x(1)],[point2y point0y(1)],'Color','c','LineWidth',2);
end
% delete(L)   % delete lines

% Get the points inside the polygon
inpoly = inpolygon(handles.featS1,handles.featS2,[point0x point0x(1)],[point0y point0y(1)]);

% Get active cell
active_cell = get(handles.popupmenu_activecell,'Value');
acnum = handles.cell_codes(active_cell);
if get(handles.radiobutton_replace,'Value')==1
    handles.CellNumbers(handles.CellNumbers==acnum) = 0;
    handles.CellNumbers(inpoly) = acnum;
elseif get(handles.radiobutton_add,'Value')==1
    handles.CellNumbers(inpoly) = acnum;
elseif get(handles.radiobutton_subtract,'Value')==1
    if ~isequal(size(handles.CellNumbers),size(inpoly))
        inpoly = inpoly';
    end
    handles.CellNumbers(handles.CellNumbers==acnum&inpoly) = 0;
end

% Clear axes
ach = allchild(handles.axes1);
delete(ach)

% Replot
replotsub(handles)

% Set active cell
popupmenu_activecell_Callback(hObject, eventdata, handles)

% Set buttondown function
bdf = 'czspikesort(''Axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function analysis_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function placefield_Callback(hObject, eventdata, handles)

% Get active cell
active_cell = get(handles.popupmenu_activecell,'Value');
acnum = handles.cell_codes(active_cell);
inx = find(handles.CellNumbers==acnum);

% Time stamps
neuron_pyr = (handles.TimeStamps(inx) - handles.TimeStamps(1)) / 1000000 + handles.corrv;

% Read position data
load([handles.inpdir '\x1_pos.mat']);
x_pos = data';
load([handles.inpdir '\x1_pos_ts.mat']);
x_pos_ts = data';
load([handles.inpdir '\y1_pos.mat']);
y_pos = data';

% Rate maps
[irhst_pyr rhst_pyr thst] = lczplace(neuron_pyr,x_pos,y_pos,x_pos_ts);
Hpyr = figure;      % plot pyramidal cell rate map
pcolor(nanpad(irhst_pyr,1))
shading flat;
grid on
colorbar
set(gca,'Layer','top','LineWidth',2,'XColor',[0.3 0.3 0.3],'YColor',...
    [0.3 0.3 0.3],'XTickLabel',{},'YTickLabel',{})

% Store rate map
eval(['handles.rhst' num2str(acnum) ' = rhst_pyr;']);
guidata(handles.figure1,handles)

% -------------------------------------------------------------------------
function [irhst,rhst2,thst] = lczplace(ncu,x_pos,y_pos,x_pos_ts)
%CZPLACE2   Place rate map.
%   [IRM RM OM] = CZPLACE2(VD,XPOS,YPOS,TS) returns place rate map (RM), 
%   interpolated place rate map (IRM, for visualization) and occupancy map
%   (OM) for discriminated unit (VD). Input arguments XPOS and YPOS are the
%   rat position coordinates, TS should contain the corresponding timestamps
%   for the position values. Bins with less than 0.25 time or 3 visits are
%   discarded.
%
%   In CZPLACE2, both the occupancy map and the rate map are smoothed by a
%   Gaussian kernel.
%
%   See also CZPLACE.

% Input argumnet check
error(nargchk(4,4,nargin))

% Position vector correction
dx = diff(x_pos);       % leave out the outliers
fdx = find(abs(dx)>8&abs([dx(2:end) 0])>8&abs(dx-[dx(2:end) 0])>16);
x_pos(fdx+1) = [];
y_pos(fdx+1) = [];
x_pos_ts(fdx+1) = [];
dy = diff(y_pos);
fdy = find(abs(dy)>8&abs([dy(2:end) 0])>8&abs(dy-[dy(2:end) 0])>16);
x_pos(fdy+1) = [];
y_pos(fdy+1) = [];
x_pos_ts(fdy+1) = [];

inxs = x_pos > 0 & y_pos > 0;   % leave out (0;0) points
x_pos2 = x_pos(inxs);
y_pos2 = y_pos(inxs);
pos_ts = x_pos_ts(inxs);

% Determine space bins
fminx = min(x_pos2);
fmaxx = max(x_pos2);
fminy = min(y_pos2);
fmaxy = max(y_pos2);

% cst = (fmaxx - fminx) / 25.5;
cst = 8;
xedge = fminx-0.0001:cst:fmaxx+cst;      % spatial bins
yedge = fminy-0.0001:cst:fmaxy+cst;
xbins = length(xedge) - 1;
ybins = length(yedge) - 1;

% Determine the position of action potentials
ncu = ncu(ncu>pos_ts(1)&ncu<pos_ts(end));     % leave action potentials before the 1st pos. timestamp
len = length(ncu);
ap_xpos = zeros(1,len);
ap_ypos = zeros(1,len);
for k = 1:len    % interpolate action potential position
    tap = ncu(k);
    ind = find(pos_ts<tap,1,'last');
    tlow = pos_ts(ind);
    xlow = x_pos2(ind);
    ylow = y_pos2(ind);
    thigh = pos_ts(ind+1);
    xhigh = x_pos2(ind+1);
    yhigh = y_pos2(ind+1);
    mpl = (tap - tlow) / (thigh - tlow);
    ap_xpos(k) = xlow + (xhigh - xlow) * mpl;
    ap_ypos(k) = ylow + (yhigh - ylow) * mpl;
end

hst = hist3([ap_xpos' ap_ypos'],'Edges',{xedge yedge})';     % spike number in space
hst = hst(1:end-1,1:end-1);     % last bin corresponds to the upper edge value

% Calculate the time spent in each space bin
thst = zeros(size(hst));        % time spent in each bin (occupancy map)
tvhst = zeros(size(hst));       % number of visits
for xb = 1:xbins
    for yb = 1:ybins
        xedlow = xedge(xb);
        xedhigh = xedge(xb+1);
        xlin = valuecrossing(pos_ts,x_pos2,xedlow,'up');
        xhin = valuecrossing(pos_ts,x_pos2,xedhigh,'down');
        xlout = valuecrossing(pos_ts,x_pos2,xedlow,'down');
        xhout = valuecrossing(pos_ts,x_pos2,xedhigh,'up');
        xin = union(xlin,xhin);
        xout = union(xlout,xhout);
        if xout(1) < xin(1)
            xin = [pos_ts(1) xin];
        end
        if xin(end) > xout(end)
            xout = [xout pos_ts(end)];
        end
        if ~isequal(length(xin),length(xout))
            error('Technical error 119.')
        end
                
        yedlow = yedge(yb);
        yedhigh = yedge(yb+1);
        ylin = valuecrossing(pos_ts,y_pos2,yedlow,'up');
        yhin = valuecrossing(pos_ts,y_pos2,yedhigh,'down');
        ylout = valuecrossing(pos_ts,y_pos2,yedlow,'down');
        yhout = valuecrossing(pos_ts,y_pos2,yedhigh,'up');
        yin = union(ylin,yhin);
        yout = union(ylout,yhout);
        if yout(1) < yin(1)
            yin = [pos_ts(1) yin];
        end
        if yin(end) > yout(end)
            yout = [yout pos_ts(end)];
        end
        if ~isequal(length(yin),length(yout))
            error('Technical error 137.')
        end
        
        for k1 = 1:length(xin)
            cxin = xin(k1);
            cxout = xout(find(xout>cxin,1,'first'));
            logix = yin<cxout&yout>cxin;
            cyins = yin(logix);
            cyouts = yout(logix);
            for k2 = 1:length(cyins)
                cyin = cyins(k2);
                cyout = cyouts(k2);
                thst(yb,xb) = thst(yb,xb) + (min(cyout,cxout) - max(cyin,cxin));
                tvhst(yb,xb) = tvhst(yb,xb) + 1;
            end
        end
    end
end
gk = gausskernel([3 3],1);
thst = smooth2_nonnan(thst,gk); % smoothing
rhst = hst ./ thst;

rhst2 = rhst .* zero2nan(double(thst>=0.25)) .* zero2nan(double(tvhst>=3));       % minimum 3 visits, 0.25 s spent in bin
% rhst2(13:15,21:23) = NaN;
rhst2 = smooth2_nonnan(rhst2,gk);
irhst = interp2_nonnan(rhst2,5);

% Plot result
figure
pcolor(rhst)
sh = findobj(allchild(gcf),'Type','surface');
set(sh,'EdgeColor','white');
c2 = [1 1 0.502; 0.9843 0.8588 0.4745; 1 0 0; 0 0.651 0; 0.502 0.502 1; 0.4157 0 0.8353];
colormap(c2)

% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles)

% Export clustering
CellNumbers = handles.CellNumbers;
DataPoints = handles.DataPoints;
TimeStamps = handles.TimeStamps;
global FEATURE_NAMES
FeatureNames = FEATURE_NAMES;
global FEATURE_VECTOR
FeatureVector = FEATURE_VECTOR;

% Save clustering variables
resdir = [handles.inpdir 'Recluster\'];     % result directory
if ~isdir(resdir)
    mkdir(resdir)
end
fn = [resdir 'clusters.mat'];       % result file
save(fn,'CellNumbers','DataPoints','TimeStamps','FeatureNames','FeatureVector')

% Export clustering quality measures
ncols = length(handles.cell_codes);
cell_codes = cell(1,ncols);
ap_num = cell(1,ncols);
isolation_dist = cell(1,ncols);
L_ratio = cell(1,ncols);
for k = 1:ncols
    currcell = handles.cell_codes(k);
    cell_codes{k} = currcell;
    ap_num{k} = length(find(handles.CellNumbers==currcell));
    [DM LrC] = isdist(k,handles);
    isolation_dist{k} = DM;
    L_ratio{k} = LrC;
end
xlsname = [resdir 'cluster_quality.xls'];   % write results to excel file
xlswrite(xlsname,cell_codes,'sheet1','A1')
xlswrite(xlsname,ap_num,'sheet1','A2')
xlswrite(xlsname,isolation_dist,'sheet1','A3')
xlswrite(xlsname,L_ratio,'sheet1','A4')
disp('Ready')

% --------------------------------------------------------------------
function extended_export_Callback(hObject, eventdata, handles)

% Export clustering
CellNumbers = handles.CellNumbers;
DataPoints = handles.DataPoints;
TimeStamps = handles.TimeStamps;
global FEATURE_NAMES
FeatureNames = FEATURE_NAMES;
global FEATURE_VECTOR
FeatureVector = FEATURE_VECTOR;

% Save clustering variables
resdir = [handles.inpdir 'Recluster_new\'];     % result directory
if ~isdir(resdir)
    mkdir(resdir)
end
fn = [resdir 'clusters.mat'];       % result file
save(fn,'CellNumbers','DataPoints','TimeStamps','FeatureNames','FeatureVector')

% Export clustering quality measures
ncols = length(handles.cell_codes);
cell_codes = cell(1,ncols);
ap_num = cell(1,ncols);
isolation_dist = cell(1,ncols);
L_ratio = cell(1,ncols);
for k = 1:ncols
    currcell = handles.cell_codes(k);
    cell_codes{k} = currcell;
    ap_num{k} = length(find(handles.CellNumbers==currcell));
    [DM LrC] = isdist(k,handles);
    isolation_dist{k} = DM;
    L_ratio{k} = LrC;
end
xlsname = [resdir 'cluster_quality.xls'];   % write results to excel file
xlswrite(xlsname,cell_codes,'sheet1','A1')
xlswrite(xlsname,ap_num,'sheet1','A2')
xlswrite(xlsname,isolation_dist,'sheet1','A3')
xlswrite(xlsname,L_ratio,'sheet1','A4')

% Save place fields
active_cell = get(handles.popupmenu_activecell,'Value');
celln = length(handles.cell_codes);
for k = 1:celln
    set(handles.popupmenu_activecell,'Value',k)
    placefield_Callback(hObject, eventdata, handles)
    handles = guidata(handles.figure1);
    fn = [resdir 'placefield' num2str(handles.cell_codes(k)) '.fig'];
    saveas(gcf,fn)
end
set(handles.popupmenu_activecell,'Value',active_cell)
acnum = handles.cell_codes(active_cell);

% Complementarity index, linear correlation (place field similarity)
for k = 1:celln
    eval(['rhst{k} = handles.rhst' num2str(handles.cell_codes(k)) ';']);
end
fn = [resdir 'placefields.mat'];       % result file
save(fn,'rhst')
pairn = size(handles.pair_codes,1);
for k = 1:pairn
    intinx = find(handles.cell_codes==handles.pair_codes(k,1));
    s_int = rhst{intinx};
    ms3_int = s_int .* (zero2nan(double(s_int>b_max_nonnan(s_int)*0.5)));
    pyrinx = find(handles.cell_codes==handles.pair_codes(k,2));
    s_pyr = rhst{pyrinx};
    ms3_pyr = s_pyr .* (zero2nan(double(s_pyr>b_max_nonnan(s_pyr)*0.5)));
    c = length(s_int(~isnan(s_int)));
    [C Cb Cc] = czcmpl2(ms3_int,ms3_pyr,c);
    R = czpfs(s_int,s_pyr);
    Rmod = czpfs_mod(s_int,s_pyr);
    fn = [resdir 'complementarity' num2str(handles.pair_codes(k,1)) '_' ...
        num2str(handles.pair_codes(k,2)) '.mat'];       % result file
    save(fn,'C','Cb','Cc','R','Rmod')
end

% Autocorrelation
for k = 1:celln
    inx = find(handles.CellNumbers==handles.cell_codes(k));
    vdisc = (handles.TimeStamps(inx) - handles.TimeStamps(1)) / 1000000;
    czacorr(vdisc)
    A_autoint = gca;
    figure
    lczacorr(vdisc)
    A_autoint2 = gca;
    
    ach1 = findobj(allchild(A_autoint),'Type','line');  % interneuron autocorrelation
    ach2 = findobj(allchild(A_autoint),'Type','hggroup');
    set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
    set(A_autoint,'XLim',[-200 200])
    axes(A_autoint)
    axis off
    fns = [resdir 'AUTO_' num2str(handles.cell_codes(k)) '.fig'];
    saveas(gcf,fns)    % save
    ach1 = findobj(allchild(A_autoint2),'Type','line');  % interneuron autocorrelation
    ach2 = findobj(allchild(A_autoint2),'Type','hggroup');
    set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
    set(A_autoint2,'XLim',[-50 50])
    axes(A_autoint2)
    axis off
    fns = [resdir 'AUTO2_' num2str(handles.cell_codes(k)) '.fig'];
    saveas(gcf,fns)    % save
end

% Normalized cross-correlation
for k = 1:pairn
    inx = find(handles.CellNumbers==handles.pair_codes(k,2));
    vdisc_x = (handles.TimeStamps(inx) - handles.TimeStamps(1)) / 1000000;
    inx = find(handles.CellNumbers==handles.pair_codes(k,1));
    vdisc_y = (handles.TimeStamps(inx) - handles.TimeStamps(1)) / 1000000;
    [H1 H2] = lczxcorr(vdisc_x',vdisc_y');   % window: +-50 ms
    figure(H1)
    A_cross = gca;
    ach1 = findobj(allchild(A_cross),'Type','line');    % crosscorrelation
    ach2 = findobj(allchild(A_cross),'Type','hggroup');
    set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
    set(A_cross,'XLim',[-50 50])
    xd = get(ach2,'XData');
    yd = get(ach2,'YData');
    axes(A_cross)
    axis off
    fns = [resdir 'CROSS_' num2str(handles.pair_codes(k,2)) '_' ...
        num2str(handles.pair_codes(k,1)) '.fig'];
    saveas(gcf,fns)    % save
    figure(H2)
    A_normcross = gca;
    ach1 = findobj(allchild(A_normcross),'Type','line');    % crosscorrelation
    ach2 = findobj(allchild(A_normcross),'Type','hggroup');
    set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
    yd = get(ach2,'YData');
    trsr = yd(52) + yd(53) + yd(54);
    set(A_normcross,'XLim',[-50 50])
    xd = get(ach2,'XData');
    yd = get(ach2,'YData');
    xdl = get(ach1(3),'XData');
    ydl = get(ach1(3),'YData');
    axes(A_normcross)
    axis off
    fns = [resdir 'NCROSS_' num2str(handles.pair_codes(k,2)) '_' ...
        num2str(handles.pair_codes(k,1)) '.fig'];
    saveas(gcf,fns)    % save
    fn = [resdir 'trsr' num2str(handles.pair_codes(k,1)) '_' ...
        num2str(handles.pair_codes(k,2)) '.mat'];       % result file
    save(fn,'trsr')
    
    [H1 H2] = lczxcorr2(vdisc_x',vdisc_y');   % window: +-5 ms
    figure(H1)
    A_smallcross = gca;
    ach1 = findobj(allchild(A_smallcross),'Type','line');    % crosscorrelation
    ach2 = findobj(allchild(A_smallcross),'Type','hggroup');
    set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
    set(A_smallcross,'XLim',[-5 5])
    xd = get(ach2,'XData');
    yd = get(ach2,'YData');
    axes(A_smallcross)
    axis off
    fns = [resdir 'SCROSS_' num2str(handles.pair_codes(k,2)) '_' ...
        num2str(handles.pair_codes(k,1)) '.fig'];
    saveas(gcf,fns)    % save
end

% Phase
for k = 1:celln
    inx = find(handles.CellNumbers==handles.cell_codes(k));
    vdisc = (handles.TimeStamps(inx) - handles.TimeStamps(1)) / 1000000 + handles.corrv;
    [H angs p_rayleigh p_rao] = lczphase(vdisc,handles);
    fns = [resdir 'PHASEHIST' num2str(handles.cell_codes(k)) '.fig'];
    saveas(H,fns)
    fn = [resdir 'PHASESTAT' num2str(handles.cell_codes(k)) '.mat'];
    save(fn,'p_rayleigh','p_rao','angs')
end

% -------------------------------------------------------------------------
function [H angs p_rayleigh p_rao] = lczphase(vdisc,handles)
%CZPHASE   Theta phase for hippocampal interneurons.
%   CZPHASE plots and saves theta phase histograms for hippocampal
%   interneurons. Theta phase is calculated via Hilbert-transform of the EEG
%   filtered between 4 and 12 Hz. First trigonometric moments are also
%   saved for all cells.
%
%   See also CZRIPPLE.

% Load
global DATAPATH
inpdir_eeg = [DATAPATH 'Czurko\czurko_EEG\'];
xlsname = [inpdir_eeg 'EEG2.xls'];
[ntz mtz] = xlsread(xlsname,'all');
o = find(strcmp(handles.filename,mtz(:,1)));
o = o(1);
fn = [inpdir_eeg 'EEG_' mtz{o,3} '_' mtz{o,1} '.mat'];
load(fn)        % load EEG
eval(['Eeg = ' mtz{o,3} ';']);
eval(['clear ' mtz{o,3}]);
eeg = Eeg.values;
sr = 1 / Eeg.interval;
eeg_start = Eeg.start;
% eeg_start = 0;
eeg_end = eeg_start + (Eeg.length - 1) * Eeg.interval;
eeg_times = eeg_start:Eeg.interval:eeg_end;
    
% Downsample EEG
eeg2 = eeg(1:5:end);
eeg_times2 = eeg_times(1:5:end);
sr2 = sr / 5;
        
% Phase calculation
[hang, hmvl, ftm, angs] = thetaphase(eeg2,vdisc,sr2,eeg_times2);
[Z,p_rayleigh,U,p_rao] = b_rao(angs);

% Plot
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm,xout] = histc(angs*180/pi,edges);   % phase histogram
nm = nm(1:end-1);
H = figure;
B = bar(cnts,nm'/length(angs));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(length(angs))];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
ts = [mtz{o,1} ' ' mtz{o,2}];
ts(ts=='_') = ' ';
title(ts)

% -------------------------------------------------------------------------
function [hang, hmvl, ftm, bahee] = thetaphase(eeg,vdisc,sr,eeg_times)

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[4 12]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 83 ms
% 3. discard cicles longer then 250 ms
fn = find(-diff(ahee)>2*pi-0.3);
sd = std(feeg);
inx = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    stim = eeg_times(fn(k):fn(k+1));
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) < 0.083 * sr) || (fn(k+1) - fn(k) > 0.25 * sr)
        inx = [inx find(vdisc>stim(1)&vdisc<stim(end))];
    end
end
vdisc(inx) = [];

% Phase angles - Hilbert
vdisc(vdisc<eeg_times(1)|vdisc>eeg_times(end)) = [];
lvd = length(vdisc)
bahee = zeros(1,lvd);
for k = 1:lvd    % interpolate
    k
    inx1 = find(eeg_times<vdisc(k),1,'last');
    inx2 = inx1 + 1;
    rto = (vdisc(k) - eeg_times(inx1)) / (eeg_times(inx2) - eeg_times(inx1));
    if rto > 1
        error('Technical error 55!')
    end
    bahee(k) = ahee(inx1) + (ahee(inx2) - ahee(inx1)) * rto;
end
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% -------------------------------------------------------------------------
function pushbutton_correction_Callback(hObject, eventdata, handles)

% Correction limit
correction_limit = str2num(get(handles.edit_correctionlimit,'String'));

% Get active cell
cellindex = get(handles.popupmenu_activecell,'Value');

% Cluster indeces
n = length(handles.CellNumbers);
inx = find(handles.CellNumbers==handles.cell_codes(cellindex));
cinx = find(handles.CellNumbers~=handles.cell_codes(cellindex));

% Feature matrix
global FEATURE_NAMES
X = [];
for k = 1:length(FEATURE_NAMES)
    featinx = find(strcmp(FEATURE_NAMES{k},handles.feature_field_names(:,1)));
    currfeat = handles.feature_field_names{featinx,2};
    eval(['currfeat_vector = handles.' currfeat ';']);
    if size(currfeat_vector,1) ~= 1
        currfeat_vector = currfeat_vector';
    end
    X = [X; currfeat_vector];
end

% Mean and covariance
XC = X(:,inx);
muC = mean(XC,2);
SC = cov(XC');
V = X - repmat(muC,1,n);
SCinv = inv(SC);

% Isolation Distance
D = mahal(X',XC');
DC = D(cinx);
sDC = sort(DC);
linx = length(inx);
if linx <= length(sDC)
    DM = sDC(length(inx));
else
    DM = 0;
end

% L-ratio
df = size(X,1);
aLC = 1 - chi2cdf(DC,df);
ninx = find(aLC>correction_limit);

% Correction
acnum = handles.cell_codes(cellindex);
handles.CellNumbers(cinx(ninx)) = acnum;
guidata(handles.figure1,handles)

% Clear axes
ach = allchild(handles.axes1);
delete(ach)

% Replot
replotsub(handles)

% Set active cell
popupmenu_activecell_Callback(hObject, eventdata, handles)

% Set buttondown function
bdf = 'czspikesort(''Axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)
guidata(handles.figure1,handles)

% -------------------------------------------------------------------------
function edit_correctionlimit_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit_correctionlimit_CreateFcn(hObject, eventdata, handles)

% Edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function pushbutton_correction2_Callback(hObject, eventdata, handles)

% Get active cell
active_cell = get(handles.popupmenu_activecell,'Value');
acnum = handles.cell_codes(active_cell);
cellinx = find(handles.CellNumbers==acnum);

% Time stamps
ts = handles.TimeStamps(cellinx);
ts = ts / 1000000;
isi = diff(ts);
tsm = find(isi<0.002);

% Correction
for k = 1:length(tsm)
    inx1 = cellinx(tsm(k));
    inx2 = cellinx(tsm(k)+1);
    if handles.CellNumbers0(inx1) ~= acnum
        handles.CellNumbers(inx1) = 0;
    end
    if handles.CellNumbers0(inx2) ~= acnum
        handles.CellNumbers(inx2) = 0;
    end
end

% Clear axes
ach = allchild(handles.axes1);
delete(ach)

% Replot
replotsub(handles)

% Set active cell
popupmenu_activecell_Callback(hObject, eventdata, handles)

% Set buttondown function
bdf = 'czspikesort(''Axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.axes1,'ButtonDownFcn',bdf)
ach = allchild(handles.axes1);
set(ach,'ButtonDownFcn',bdf)
guidata(handles.figure1,handles)