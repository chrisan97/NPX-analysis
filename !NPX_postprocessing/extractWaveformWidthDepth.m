function waveformWidth = extractWaveformWidthDepth(sp, goodClusters, myKsDir)
%
% This code is used when you want to extract the waveform width and Depth for a given
% cluster.
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/29/2024
%
% ------
% Inputs
% ------
% sp: spike structure holding array of information. Usually generated from
% loadKsDir
% goodClusters: array holding cluster numbers that are labeled 'good'
% visStimTimes: visual stimulus timing in array
% myKsDir: directory to the kilosort folder
%
% ------
% Outputs
% ------
% Generates an array holding waveform width for each good cluster
% Generates an array holding waveform depth for each good cluster
%

%% Iterate throough good clusters
waveformWidth = zeros(2,length(goodClusters));
waveformDepth = zeros(2,length(goodClusters));
for i = 1:length(goodClusters) % for all the good clusters in my recording

    % First, setup the classical structure for extracting waveforms
    myKsDir = cd;
    bindataDir = fullfile(cd, '..');
    gwfparams.bindataDir = bindataDir;
    gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
    apD = dir(fullfile(bindataDir, '*ap*.bin')); % AP band file from spikeGLX specifically
    gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
    gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
    gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
    gwfparams.wfWin = [-50 51];              % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf = 500;                    % Number of waveforms per unit to pull out
    gwfparams.spikeTimesSamples = sp.ss(sp.clu == goodClusters(i));
    gwfparams.spikeTimesSeconds = sp.st(sp.clu == goodClusters(i));
    gwfparams.spikeClusters = sp.clu(sp.clu==goodClusters(i));

    stLength = numel(gwfparams.spikeTimesSamples);  
    if stLength > gwfparams.nWf % If this cluster has more than 500 spikes
        thisWaveformNumber = randperm(stLength,gwfparams.nWf); % we select 500 spikes
    else
        thisWaveformNumber = randperm(stLength); % otherwise we just shuffle all the waveforms
    end

    wf = getThisWaveform(gwfparams,thisWaveformNumber);
    wfExtracted = wf.waveFormsMean(wf.maxChannel,:);
    [minVal troughIdx] = min(wfExtracted); 
    [maxVal peakIdx] = max(wfExtracted(troughIdx:end));
    wfWidth = (peakIdx-1) * 1/30000 * 1000; % convert into ms

    %% save values for width
    waveformWidth(1,i) = goodClusters(i);
    waveformWidth(2,i) = wfWidth;

    %% save values for depth
    xCoord = sp.xcoords(wf.maxChannel);
    yCoord = sp.ycoords(wf.maxChannel);
    waveformDepth(1,i) = goodClusters(i)+1;
    waveformDepth(2,i) = xCoord;
    waveformDepth(3,i) = yCoord;
    waveformDepth_Table = array2table(waveformDepth');
    waveformDepth_Table.Properties.VariableNames = {'Cluster','xCoord','yCoord'};
end

%% Save the outputs to a folder for width
fname = fullfile(myKsDir,'analysis_variables','waveformWidth');
save(fname,'waveformWidth')

%% Save the outputs toa folder for depth
xCoordsRange = [max(sp.xcoords) min(sp.xcoords)];
yCoordsRange = [max(sp.ycoords) min(sp.ycoords)];
fname = fullfile(myKsDir,'analysis_variables','waveformDepth');
save(fname,'waveformDepth_Table','xCoordsRange','yCoordsRange')

end