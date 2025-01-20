function sAmp = getSpikeAmplitudes(gsaparams)
% function sAmp = getSpikeAmplitudes(gwfparams)
%
% Extracts individual spike waveforms from the raw datafile and grabs its
% amplitude across time
%
% Modified from getWaveForms() by Chris (Seong Yeol) An
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimesSamples =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
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
fileName = fullfile(gsaparams.bindataDir,gsaparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gsaparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gsaparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gsaparams.wfWin(1):gsaparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gsaparams.dataType, [gsaparams.nCh nSamp], 'x'});
chMap = readNPY(fullfile(gsaparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

% Read spike time-centered waveforms
unitIDs = unique(gsaparams.spikeClusters);
numUnits = size(unitIDs,1);

allAmplitude = nan(numUnits,gsaparams.lenOfSpikeTrain-1); % Here I am taking the last spike out (to make sure waveform doesnt get cut off)

for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gsaparams.spikeTimesSamples(gsaparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    for i = 1:curUnitnSpikes-1
        tmpWf = mmf.Data.x(gsaparams.maxChannel, curSpikeTimes(i)+gsaparams.wfWin(1):curSpikeTimes(i)+gsaparams.wfWin(2));
        [minVal minCol] = min(tmpWf);
        [maxVal maxCol] = max(tmpWf(minCol:end));         
        allAmplitude(i) = maxVal - minVal;       % Find amplitude of the spike   
    end
end

% Package in wf struct
sAmp.unitIDs = unitIDs;
sAmp.amplitude = allAmplitude;


end
