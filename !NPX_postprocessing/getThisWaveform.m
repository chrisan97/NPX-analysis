function wf = getThisWaveform(gwfparams,thisWaveformNumber)
% function wf = getThisWaveform(gwfparams)
%
% Grab the waveform at the given timepoint
%
% Contributed by C. Schoonover and A. Fink
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.bindataDir,gwfparams.fileName);  
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

%% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
waveForms = nan(gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(nChInMap,wfNSamples);

curUnitID = unitIDs;
curSpikeTimes = gwfparams.spikeTimesSamples(gwfparams.spikeClusters==curUnitID);
curSpikeTimes_noEnd = curSpikeTimes(1:end-2);
curUnitnSpikes = size(curSpikeTimes_noEnd,1);
spikeTimeKeeps = curSpikeTimes(thisWaveformNumber);
for curSpikeTime = 1:length(spikeTimeKeeps)
    tmpWf = mmf.Data.x(1:nChInMap,spikeTimeKeeps(curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curSpikeTime)+gwfparams.wfWin(end));
    waveForms(curSpikeTime,:,:) = tmpWf(chMap,:);
end
waveFormsMean = squeeze(nanmean(waveForms(:,:,:),1));
squeezedWf = (waveFormsMean * 0.5 / 8192) / 80 * 1000000;     

%% Find max channel for each waveform
wf.maxChannel = [];
totalWaveform = [];
for j = 1:nChInMap
    [minVal minCol] = min(squeezedWf(j,1:end)); 
    [maxVal maxCol] = max(squeezedWf(j,minCol:end));
    totalWaveform(j,1) = maxVal - minVal;       % Find amplitude of the spike
    totalWaveform(j,2) = (maxCol - 1) * 1/30000 * 1000;      % Because maxCol index starts from the index of minCol
end
[maxAmplitudeWaveform maxChannel] = max(totalWaveform);     % Take the waveform with maximum amplitude and grab the maximum channel 
wf.maxChannel = maxChannel(1);

%% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;
end
