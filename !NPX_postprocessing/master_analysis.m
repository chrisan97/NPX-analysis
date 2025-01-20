%% Analysis master file
% Run these functions in order to generate given plots. 
%
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/21/2024
%
%
%% Set up paths to import the data
format longg 
clear all
clc

myKsDir = cd; % cd to the kilosort output folder
Ksparams.excludeNoise = false; % we make sure all the spike times are imported. 
sp = loadKSdir(myKsDir,Ksparams); % we load the structure in.
spiketimes = readmatrix('spike_times.txt');  % this is converted using the accurate spike rate time
sp.st = spiketimes; 
sp.ss = readNPY('spike_times.npy'); % in samples

goodClu_idx = find(sp.cgs == 2); % find index of good clusters
goodClusters = sp.cids(goodClu_idx); % this is -1 from the MATLAB indexing

%% Load visual stimulation timestamps
bindataDir = fullfile(cd, '..','*_1500_TPrime.txt');
visStimFiles = dir(bindataDir);
visStimTimes = readmatrix(fullfile(visStimFiles.folder,visStimFiles.name));

%% Load the order of visual stimulation
eventDir= fullfile(cd, '..','*_visStimOrder.txt');
eventDirFiles = dir(eventDir);
eventTable = readtable(fullfile(eventDirFiles.folder,eventDirFiles.name),'Delimiter',{':','-'});
eventTable.Var3(strcmp(eventTable.Var3,'')) = {NaN}; % these lines are
%sometimes needed when the textfile is differently formatted, unsure why..
eventTable.Var3 = str2double(eventTable.Var3);
 
RPiVisAllMap = eventTable(contains(eventTable.Var1,'Sine map'),:);
RPiVisAllMap.StimTimes = visStimTimes(1:height(RPiVisAllMap));

RPiVisSineMap = RPiVisAllMap(matches(RPiVisAllMap.Var1,'Sine map'),:); % these tables may be empty depending on the stimulation protocol
RPiVisGratingMap = RPiVisAllMap(matches(RPiVisAllMap.Var1,'Grating map'),:); % these tables may be empty depending on the stimulation protocol
RPiVisSquareMap = RPiVisAllMap(matches(RPiVisAllMap.Var1,'Square map'),:); % these tables may be empty depending on the stimulation protocol

%% Load "blank map" which is 500ms stimulation 
PowerLevelMap = eventTable(contains(eventTable.Var1,'Blank map'),:);

% We extrac optoStimTimes. The user need to provide indices to
% include/exclude from optoTimes. The below times are customized
% 1:200 optoTag
% 201:260 optoSine
% 261:284 10.5mW
% 285:308 5mW
% 309:332 8mW
% 333:356 3mW
% 357:380 15mW
% 381:580 optoTag
mw_10_idx = 261:284;
mw_10_optoTimes = optoDuration.optoTimes(mw_10_idx);

mw_5_idx = 285:308;
mw_5_optoTimes = optoDuration.optoTimes(mw_5_idx);

mw_8_idx = 309:332;
mw_8_optoTimes = optoDuration.optoTimes(mw_8_idx);

mw_3_idx = 333:356;
mw_3_optoTimes = optoDuration.optoTimes(mw_3_idx);

mw_15_idx = 357:380;
mw_15_optoTimes = optoDuration.optoTimes(mw_15_idx);

% combine optotimes
mw_optoTimes = [mw_10_optoTimes; mw_5_optoTimes; mw_8_optoTimes; mw_3_optoTimes; mw_15_optoTimes]; % this only works when they are same length
mw_assignment = [repmat(10.5,length(mw_10_optoTimes),1); repmat(5,length(mw_5_optoTimes),1); repmat(8,length(mw_8_optoTimes),1);...
    repmat(3,length(mw_3_optoTimes),1); repmat(15,length(mw_15_optoTimes),1)];

mw_Table = table(mw_assignment,mw_optoTimes,'VariableNames',{'Power','OptoTimes'});

assert(length(mw_optoTimes) == sum(PowerLevelMap.Var3 == 1))



%% Load optogenetics times
optoRiseDir = fullfile(cd, '..', '*xd_3_7_0_TPrime.txt');
optoDirFiles = dir(optoRiseDir);
optoRiseTable = readtable(fullfile(optoDirFiles.folder,optoDirFiles.name));

optoFallDir = fullfile(cd, '..', '*xid_3_7_0_TPrime.txt');
optoDirFiles = dir(optoFallDir);
optoFallTable = readtable(fullfile(optoDirFiles.folder,optoDirFiles.name));

assert(height(optoRiseTable) == height(optoFallTable)) % The number of rising edges and falling edges should be the same.

% Calculate Opto diff
optoDuration = optoFallTable - optoRiseTable;

% Now I want to calculate the stimulus type used for optogenetics
% stimulation
StimType = {};
for i = 1:height(optoDuration)
    if optoDuration.Var1(i) >= 0.45 & optoDuration.Var1(i) <= 0.55 % If around 500ms, we assume that it is 500ms. 
        StimType{end+1} = '0.500';
    elseif optoDuration.Var1(i) >= 0.0028 & optoDuration.Var1(i) <= 0.0037 % if around 3.3ms, we assume that is is 3ms
        StimType{end+1} = '0.0033';
    elseif optoDuration.Var1(i) >= 0.0018 &  optoDuration.Var1(i) <= 0.00279
        StimType{end+1} = '0.002';
    else
        StimType{end+1} = 'unclassified';
    end
end

optoDuration = addvars(optoDuration,StimType');
optoTimes = optoRiseTable.Var1;
optoDuration = addvars(optoDuration,optoTimes);
optoConds = unique(optoDuration.Var2); % These are list of optogenetics stimulation conditions


%% Functions to generate individual subplots

% This will generate scatters and plots considering orientation (saved in folder Orientation)
plotVisStimScatter(sp, goodClusters, myKsDir, RPiVisAllMap)
close all;

% This will generate scatters and plots considering orientation (with
% optogenetics) 
% !! run this if you have visual stimulus + 500ms opotogenetics (saved in
% folder OptoOrientation)
plotVisStimOptoScatter(sp, goodClusters, visStimTimes, myKsDir, RPiVisAllMap)
close all;

%
plotOptoScatter(sp, goodClusters, myKsDir, RPiVisAllMap, optoDuration)
close all;

% Plot OSI things
calculateOSI(goodClusters, myKsDir);

% This compares firing rate across all conditions
calculateFRComparison(sp, goodClusters, myKsDir, RPiVisAllMap)
close all

% Calculate Waveform width
waveformWidth = extractWaveformWidthDepth(sp, goodClusters, myKsDir);

calculateOptoModIndex(sp, goodClusters, myKsDir, RPiVisAllMap, optoDuration)