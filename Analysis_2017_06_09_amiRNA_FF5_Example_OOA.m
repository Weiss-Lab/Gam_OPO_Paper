%% Imports
% close all 
%#ok<*AGROW>
warning('off', 'MATLAB:xlswrite:AddSheet')
clc

% General experiment descriptors
expName = 'amiRNA_FF5_Example_OOA';
expDate = '2017_06_09';
expCyto = 'SBC-Fortessa';
expFolder = pwd;
controlsFolder = [pwd, filesep, 'Controls', filesep];

% Load all data into path, including back end code (since it is in this folder, otherwise add separately)
addpath(genpath(expFolder));

% Create analysis namespaces
flow = FlowAnalysis();
plt = Plotting();
trf = Transforms();
comp = Compensation();


fprintf(1, 'Finished adding paths\n');
saveFileName = strcat('Analysis_', expDate, '_', expName);
fprintf(1, 'Analysis save file name:\t\t%s\n', saveFileName);



%% Prepare Experiment Details

% Names of channels
channels = {'Pacific_Blue_A', 'FITC_A', 'PE_Texas_Red_A'};
BLU = channels{1};	% EBFP2
YEL = channels{2};	% EYFP
RED = channels{3};	% mKate2

% Define experiment details
expDetails = struct( ...
	'date', expDate, ...
	'name', expName, ...
	'folder', expFolder, ...
	'cytometer', expCyto);

% Define rainbow calibration beads 
beads = struct( ...
	'filename', 's_JG_beads_211.fcs', ...
	'type', 'RCP-30-5A', ...
	'lot', 'AH01', ...
	'date', expDate, ...
	'cytometer', expCyto);

% Sample filenames
sampleFnames = {'s_FT_JG412_JG365_197.fcs'};

% Control filenames
wtFname = {'s_FT_NoDNA_209.fcs'};
scFnames = {'s_FT_EBFP2_207.fcs', 's_FT_EYFP_208.fcs', 's_FT_mKate2_206.fcs'};


fprintf(1, 'Finished setting experiment details\n');



%% Initialize and Process Data

	
% Load data
flowData = FlowData(sampleFnames, channels, expDetails);

% Add controls
flowData.addControls(controlsFolder, wtFname, scFnames);

% Gate
flowData.gate();

% Compensate
flowData.compensate('raw', 'P3', struct( ...
	'plotsOn', true, 'recompute', true, 'plotLin', true));


fprintf(1, '---> Finished processing!\n\n');



%% Scatterplots

% Figure setup
numPoints = 10000;
cmap = parula(100);
figScatter = figure();
ax = gca();

% Data slicing info
dataType = 'comp';
gate = 'P3';
units = ' (AFU)';

% Slice data
outData = flowData.slice(1, struct( ...
		'channels', {channels}, ...
		'dataType', dataType, ...
		'gate',		gate));

% Logicle-transform and extract data
outDataBX = trf.lin2logicle(outData);
valid = (sum(outDataBX > 0, 2) == size(outDataBX, 2)); % Ignore neg data
bluData = outDataBX(valid, 1);
yelData = outDataBX(valid, 2);
redData = outDataBX(valid, 3);

% Subsample fewer points
ss = flow.subSample(numel(bluData), numPoints);

% Plot points in 3D scatter colored by mKate2 output level
[colors, sortIdx] = plt.getColors(redData, cmap, struct('min', 0, 'max', 4.5));
scatter3(ax, bluData(sortIdx), yelData(sortIdx), redData(sortIdx), 8, colors, 'filled')

% Apply biexponential axes tick marks/labels
plt.biexpAxes(ax, true, true, true);
ax.FontSize = 12;

% Title/axes labels
title(ax, 'miR-FF5', 'fontsize', 16)
xlabel(ax, ['EBFP2-Marker', units])
ylabel(ax, ['EYFP-miRNA', units])
zlabel(ax, ['mKate2-Ouptut', units])


fprintf(1, 'Finished plotting figures\n');



%% Bin (Full 2D)

% Define bins
edgesBlu = trf.logicle2lin(linspace(0.1, 4, 17));
edgesYel = trf.logicle2lin(linspace(0.5, 4, 9));
binInputs = struct(BLU, edgesBlu, YEL, edgesYel);

% Bin
flowData.bin(binInputs, 'comp');


fprintf(1, 'Finished binning all data\n');



%% Bin Highlighting

units = ' (AFU)';
dataType = 'comp';
gate = 'P3';
figBinColors = figure();

% Define bin colors
cmapC = cool(flowData.binSizes(1));

cmapRY = [
	150 150 150
    214 223 35
    251 175 63
    247 147 29
    240 90 40
    239 64 54
    190 30 45
    150 150 150]./255;

cmapP = parula(flowData.numBins);

% --- EBFP2 Bins --- %

ax = subplot(1, 3, 1); hold(ax, 'on');
	
% First plot all data in grey, then highlight the desired bins w/ colors
nonOutData = flowData.slice(1, struct( ...
		'channels', {{BLU, YEL}}, ...
		'dataType', dataType, ...
		'gate',		gate));

scatter(trf.lin2logicle(nonOutData(:, 2)), ...
		trf.lin2logicle(nonOutData(:, 1)), ...
		5, [0.7, 0.7, 0.7], 'filled')

for b = 1:flowData.binSizes(1)
	
	% Define bins to extract
	bins = b * ones(flowData.binSizes(2), 2);
	bins(:, 2) = 1:flowData.binSizes(2);
	
	outData = flowData.slice(1, struct( ...
			'channels', {{BLU, YEL}}, ...
			'dataType', dataType, ...
			'gate',		gate, ...
			'bins',		bins));
	
	scatter(trf.lin2logicle(outData(:, 2)), ...
			trf.lin2logicle(outData(:, 1)), ...
			5, cmapC(b, :), 'filled')
end
plt.biexpAxes(ax, true, true);
title(ax, 'Demonstration of EBFP2 Binning', 'fontsize', 16)
xlabel(ax, ['EBFP2-Marker', units])
ylabel(ax, ['EYFP-miRNA', units])
ax.FontSize = 12;


% --- EYFP Bins --- %

ax = subplot(1, 3, 2); hold(ax, 'on');
	
% First plot all data in grey, then highlight the desired bins w/ colors
nonOutData = flowData.slice(1, struct( ...
		'channels', {{BLU, YEL}}, ...
		'dataType', dataType, ...
		'gate',		gate));

scatter(trf.lin2logicle(nonOutData(:, 2)), ...
		trf.lin2logicle(nonOutData(:, 1)), ...
		5, [0.7, 0.7, 0.7], 'filled')

for b = 1:flowData.binSizes(2)
	
	% Define bins to extract
	bins = b * ones(flowData.binSizes(1), 2);
	bins(:, 1) = 1:flowData.binSizes(1);
	
	outData = flowData.slice(1, struct( ...
			'channels', {{BLU, YEL}}, ...
			'dataType', dataType, ...
			'gate',		gate, ...
			'bins',		bins));
	
	scatter(trf.lin2logicle(outData(:, 2)), ...
			trf.lin2logicle(outData(:, 1)), ...
			5, cmapRY(b, :), 'filled')
end
plt.biexpAxes(ax, true, true);
title(ax, 'Demonstration of EYFP Binning', 'fontsize', 16)
xlabel(ax, ['EBFP2-Marker', units])
ylabel(ax, ['EYFP-miRNA', units])
ax.FontSize = 12;


% --- 2D Bins --- %

ax = subplot(1, 3, 3); hold(ax, 'on');
	
% First plot all data in grey, then highlight the desired bins w/ colors
nonOutData = flowData.slice(1, struct( ...
		'channels', {{BLU, YEL}}, ...
		'dataType', dataType, ...
		'gate',		gate));

scatter(trf.lin2logicle(nonOutData(:, 2)), ...
		trf.lin2logicle(nonOutData(:, 1)), ...
		5, [0.7, 0.7, 0.7], 'filled')

randC = randperm(flowData.numBins);
for b = 1:flowData.numBins
		
	outData = flowData.slice(1, struct( ...
			'channels', {{BLU, YEL}}, ...
			'dataType', dataType, ...
			'gate',		gate, ...
			'bins',		b));
	
	scatter(trf.lin2logicle(outData(:, 2)), ...
			trf.lin2logicle(outData(:, 1)), ...
			5, cmapP(randC(b), :), 'filled')
end
plt.biexpAxes(ax, true, true);
title(ax, 'Demonstration of 2D Binning', 'fontsize', 16)
xlabel(ax, ['EBFP2-Marker', units])
ylabel(ax, ['EYFP-miRNA', units])
ax.FontSize = 12;


fprintf(1, 'Finished plotting figures\n')



%% Extract Bin Data

binData = struct();
stat = 'p50'; % Median (50th %ile)

% Extract bin data
for ch = channels
	
	% Slice parameters to extract bin statistics
	sliceParams = struct( ...
		'channels', ch{:}, ...
		'dataType', 'comp', ...
		'gate', 'P3', ...
		'bins', 'all');
	
	binData.(ch{:}) = flowData.computeBinStats( ...
			1, sliceParams, stat).(stat);
		
	% Compute normalized bin data to the low-EYFP bins (~no miR-FF5)
	binDataNorm.(ch{:}) = binData.(ch{:}) ./ binData.(ch{:})(:, 1);
end

save binData.mat binData binDataNorm


fprintf(1, 'Finished extracting bin data\n')



%% Plot 2D Bin Data Heatmaps

% load('binData.mat'); % If re-loading

% Axes properties (can be adjusted)
axProperties = struct();

% Set edges and labels for the plot
edges = {trf.lin2logicle(edgesBlu), ...
		 trf.lin2logicle(edgesYel)};
labels = {['EBFP2-Marker', units], ...
		['EYFP-miRNA', units], ...
		['mKate2-Output', units]};

% Heatmap options
options = struct( ...
	'biexp', {{'X', 'Y', 'C'}}, ...
	'min', 0, 'max', 4.5);	% Min and max specify colorbar limits

figBinHmap = plt.binHeatmap(trf.lin2logicle(binData.(RED)), ...
	edges, labels, parula(100), axProperties, options);


fprintf(1, 'Finished plotting figures\n');



%% Plot Norm 2D Bin Data Heatmap

% load('binData.mat'); % If re-loading

% Axes properties (can be adjusted)
axProperties = struct();

% Set edges and labels for the plot
edges = {trf.lin2logicle(edgesBlu), ...
		 trf.lin2logicle(edgesYel)};
labels = {['EBFP2-Marker', units], ...
		['EYFP-miRNA', units], ...
		['mKate2 Fold-\Delta', units]};

% Heatmap options
options = struct( ...
	'biexp', {{'X', 'Y'}}, ...
	'min', 0, 'max', 2.0);	% Min and max specify colorbar limits

figBinHmapNorm = plt.binHeatmap(binDataNorm.(RED), ...
	edges, labels, ColorMap('redblue').getColormap(100), axProperties, options);


fprintf(1, 'Finished plotting figures\n');



%% Plot 2D Bin Data Surface Plots

% load('binData.mat'); % If re-loading

% Figure setup
cmap = parula(100);
figBinSurf = figure();
ax = gca();

% Data slicing info
dataType = 'comp';
gate = 'P3';
units = ' (AFU)';

% Plot points in 3D scatter colored by mKate2 output level
colors = plt.getColors(trf.lin2logicle(binData.(RED)), cmap, struct('min', 0, 'max', 4.5));
surface(ax, trf.lin2logicle(binData.(BLU)), ...
		trf.lin2logicle(binData.(YEL)), ...
		trf.lin2logicle(binData.(RED)), colors)

% Apply biexponential axes tick marks/labels
plt.biexpAxes(ax, true, true, true);
ax.FontSize = 12;

% Set camera view angle to the same as for the scatterplot
ax.CameraViewAngle = 11.1414;
ax.CameraPosition = [-17.6377, -17.1149, 27.4343];

% Title/axes labels
title(ax, 'miR-FF5', 'fontsize', 16)
xlabel(ax, ['EBFP2-Marker', units])
ylabel(ax, ['EYFP-miRNA', units])
zlabel(ax, ['mKate2-Ouptut', units])


fprintf(1, 'Finished plotting figures\n');



%% Transfer Curves (Fluor)

figBinsTCs = figure();
ax = gca(); hold(ax, 'on');
units = ' (AFU)';

for b = 1:flowData.binSizes(2)
	
	plot(ax, trf.lin2logicle(binData.(BLU)(:, b)), ...
		trf.lin2logicle(binData.(RED)(:, b)), ...
		'.--', 'color', cmapRY(b, :), 'markersize', 30, 'linewidth', 2)
	
end

% Apply biexponential axes tick marks/labels
plt.biexpAxes(ax, true, true);
ax.FontSize = 12;

% Title/axes labels
title(ax, 'miR-FF5', 'fontsize', 16)
xlabel(ax, ['EBFP2-Marker', units])
ylabel(ax, ['mKate2-Ouptut', units])

legend(ax, {'Bin 1', 'Bin 2', 'Bin 3', 'Bin 4', 'Bin 5', 'Bin 6', 'Bin 7', 'Bin 8'});

fprintf(1, 'Finished plotting figures\n');







