function calculateOptoModIndex(sp, goodClusters, myKsDir, RPiVisAllMap, optoDuration)
%
% This plot is used when we want to align the spikes to the optogenetics
% pulse, while just ignoring other cases of visual presentation.
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/31/2024
%
% ------
% Inputs
% ------
% sp: spike structure holding array of information. Usually generated from
% loadKsDir
% goodClusters: array holding cluster numbers that are labeled 'good'
% visStimTimes: visual stimulus timing in array
% myKsDir: directory to the kilosort folder
% RPiVisAllMap: a table holding visual stimulation times, type of vis stim, etc 
% optoDuration: a table holding the length of optogenetic time, the
% "classification of optogenetic type", and the timestamp.
%
% ------
% Outputs
% ------
% Generates a plot for each unit.
%

%% There are two ways of calculating optoModulationIdx
% 1. We take the first 500ms of the visual stimulus as baseline and the
% next 500ms (with opto) as active
%
% 2. Since there are moving gratings, this does not guarantee that the
% baseline firing at these two epochs are equivalent. Thus, we take
% baseline as non-opto trials 500ms ~ 1000ms and active as opto trials
% 500ms ~ 1000ms.
% 

% We first want to find all times when optogenetic stimulation is presented
% concurrently with the sine grating. RPiVisAllMap.Var3 = 1 will represent
% this.
concurrentOptoVisStimIdx = (RPiVisAllMap.Var3 == 1);
concurrentOptoVisStimTimes = RPiVisAllMap.StimTimes(concurrentOptoVisStimIdx);
nonOptoVisStimIdx = (RPiVisAllMap.Var3 == 0);
nonOptoVisStimTimes = RPiVisAllMap.StimTimes(nonOptoVisStimIdx);

%% Initialize some important variables
optoConds = unique(optoDuration.Var2); % These are list of optogenetics stimulation conditions
colors = {'k','k','k'};

%% Iterate through each cluster
optoModulationMatrix = zeros(length(goodClusters),4); % Firing rate from 500-1000ms

for i = 1:length(goodClusters) % for all the good clusters in my recording
    unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds, take spikes from the good cluster
    
    % (method 1) To calculate optomodulation index, we calculate
    % Active - baseline/Active + baseline
    % We use active as optogenetics onset from 500~1000
    % for baseline, we use 0~500ms in the same trial
    tempBaseline = [];
    tempActive = [];
    tempOtherBaseline = [];
    % Calculate baselina and active FR
    for stim = 1:length(concurrentOptoVisStimTimes)
        % First, I want to find the time of the optogenetics stimulation in
        % this trial
        thisTrialOptoTime = optoDuration.optoTimes(find(optoDuration.optoTimes >= concurrentOptoVisStimTimes(stim),1));       
        alignedThisTrialOptoTime = thisTrialOptoTime - concurrentOptoVisStimTimes(stim); % aligned to the timeframe of stim
        assert(alignedThisTrialOptoTime <= 0.7); % make sure this is the correct one
       
        % Then, align the spikes to the onset of the visual stimulus
        alignSpikes = unitSpikes - concurrentOptoVisStimTimes(stim); 
        baselineSpikes = alignSpikes(alignSpikes >= 0 & alignSpikes < alignedThisTrialOptoTime); % holds spiketimes
        activeSpikes = alignSpikes(alignSpikes >= alignedThisTrialOptoTime & alignSpikes < alignedThisTrialOptoTime + 0.5); % holds spiketimes
    
        % Lastly, we convert the spikes into firing rates
        baselineFR = length(baselineSpikes)/(alignedThisTrialOptoTime); % in Hz
        activeFR = length(activeSpikes)/0.500; % in Hz  

        tempBaseline = [tempBaseline baselineFR];
        tempActive = [tempActive activeFR];
    end   
    % Calculate baseline firing rate for other interleaved non-opto trials
    for stim = 1:length(nonOptoVisStimTimes)
        alignSpikes = unitSpikes - nonOptoVisStimTimes(stim);
        baselineSpikes = alignSpikes(alignSpikes >= 0.5 & alignSpikes < 1.0);
        baselineFR = length(baselineSpikes)/0.500; % in Hz

        tempOtherBaseline = [tempOtherBaseline baselineFR];
    end
    % save into overall matrix
    optoModulationMatrix(i,1) = goodClusters(i) + 1;
    optoModulationMatrix(i,2) = sum(tempBaseline)/length(tempBaseline); % baseline
    optoModulationMatrix(i,3) = sum(tempOtherBaseline)/length(tempOtherBaseline); % other baseline
    optoModulationMatrix(i,4) = sum(tempActive)/length(tempActive); % active
end

% Convert the matrix into table
columnNames = {'Cluster','Baseline','otherBaseline','Active'};
optoModulation_Table = array2table(optoModulationMatrix,'VariableNames',columnNames);
% Calculate OMI
optoModulationIdx = (optoModulation_Table.Active - optoModulation_Table.Baseline)./(optoModulation_Table.Active + optoModulation_Table.Baseline);
optoModulationIdxOther = (optoModulation_Table.Active - optoModulation_Table.otherBaseline)./(optoModulation_Table.Active + optoModulation_Table.otherBaseline);
optoModulation_Table.OMI = optoModulationIdx;
optoModulation_Table.otherOMI = optoModulationIdxOther;

% Save the table
fname = fullfile(myKsDir,'analysis_variables','optoModulationTable');
save(fname,'optoModulation_Table')


end