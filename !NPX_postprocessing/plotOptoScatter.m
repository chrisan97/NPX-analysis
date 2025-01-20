function plotOptoScatter(sp, goodClusters, myKsDir, RPiVisAllMap, optoDuration)
%
% This plot is used when we want to align the spikes to the optogenetics
% pulse, while just ignoring other cases of visual presentation.
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/18/2024
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

%% Initialize some important variables
optoConds = unique(optoDuration.Var2); % These are list of optogenetics stimulation conditions
colors = {'k','k','k'};
binSize = 0.020;
edges = -0.3:binSize:2.0;
edgesOrientation = 0:binSize:1.5;
edgesPartial = 0.5:binSize:1.5;
numOfbins = length(edges)-1; % calculate the number of bins

%% Iterate through each cluster
for i = 1:length(goodClusters) % for all the good clusters in my recording
    unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds, take spikes from the good cluster
    numOfSpikesPSTH = zeros(length(optoConds),numOfbins); % initialize the array to hold PSTHs

    unitSpikeOptoScatter = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);
    % First plot for all opto-conditions
    for cond = 1:length(optoConds)
        rowIdx = mod(cond-1,length(optoConds)) + 1;
        colIdx = ceil(cond/length(optoConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;
        subplot(length(optoConds),4,subplotIdx)
    
        y = 0; % for plotting
        numOfStims = sum(matches(optoDuration.Var2,optoConds(cond))); % how many instances of the stimulation in this condition
        thisOptoStimTimes = optoDuration.optoTimes(matches(optoDuration.Var2,optoConds(cond)));
        assert(length(thisOptoStimTimes) == numOfStims) % just make sure extraction was correct
        stimLength = str2num(cell2mat(optoConds(cond)));

        nonOptoWaveformNumber = []; % This holds the last spike before optogenetics pulse
        firstOptoWaveformNumber = []; % These are for comparing optogenetics waveforms. This holds the first waveform after stim
        firstOptoSpikeLatency = []; % Latency of the first spike
        otherOptoWaveformNumber = []; % This holds the second waveform
        for stim = 1:numOfStims % for each stimulation in in the given condition
            alignSpikes = unitSpikes - thisOptoStimTimes(stim); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3);
            N = histcounts(indivTrial, edges);
            numOfSpikesPSTH(cond,:) = numOfSpikesPSTH(cond,:) + N;
            % Here we are looking for any spike that falls in this window
            % of 7ms
            temp_firstSpike = find(alignSpikes >= 0 & alignSpikes <= 0.007,1); 
            if isempty(temp_firstSpike)
                firstOptoWaveformNumber = [firstOptoWaveformNumber NaN];
                firstOptoSpikeLatency = [firstOptoSpikeLatency NaN];
            else
                firstOptoWaveformNumber = [firstOptoWaveformNumber temp_firstSpike];
                firstOptoSpikeLatency = [firstOptoSpikeLatency alignSpikes(temp_firstSpike)];
                scatter(alignSpikes(temp_firstSpike),y,20,'r','filled')
                hold on;
            end
            % Now we want to find any spike that comes a little later in
            % the stimulation. Assumption: first spike could have some
            % weird artifacts? So let's just find any spike between 50ms
            % and 100ms
            temp_otherSpike = find(alignSpikes >= 0.05 & alignSpikes <= 0.1,1);
            if isempty(temp_otherSpike)
                otherOptoWaveformNumber = [otherOptoWaveformNumber NaN];
            else
                otherOptoWaveformNumber = [otherOptoWaveformNumber temp_otherSpike];
            end
            % Now I want to find a spike that is outside of the stimulation
            % times. Ideally, this is a spike that is right before the
            % optogenetics stimulation. If this is not the case then we can
            % randomly draw from non-opto spikes.
            temp_lastSpike = find(alignSpikes < 0 & alignSpikes >= -0.1,1,'last'); % Here we are looking for spikes 100ms before the onset of stimulation
            if isempty(temp_lastSpike)
                nonOptoWaveformNumber = [nonOptoWaveformNumber NaN];
            else
                nonOptoWaveformNumber = [nonOptoWaveformNumber temp_lastSpike];
            end
        end
        % Plot other aspects of the scatter (line at 0, etc)
        plot(zeros(1, numOfStims), 1:numOfStims, 'k--', 'LineWidth', 1.0)
        hold on;
        fill([0, stimLength, stimLength, 0], [0, 0, numOfStims, numOfStims], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
        xlim([0-stimLength, stimLength * 3]);
        ylim([0, numOfStims]);
        xlabel('Time (s)');
        ylabel('Optogenetics pulse');
        title(append('Response to ',optoConds{cond}));
        ax = gca;
        ax.XAxis.Exponent = 0;
        ax.FontSize = 10;

        % Some important measures of opto-tagging
        firstSpikeReliability = sum(~isnan(firstOptoWaveformNumber))/numOfStims;
        firstSpikeLatency = mean(firstOptoSpikeLatency,'omitnan');
        firstSpikeJitter = std(firstOptoSpikeLatency,'omitmissing');
        xPos = stimLength * 2.8;
        yPos = numOfStims * 0.95;
        text(xPos, yPos, sprintf('Reliability: %.2f\nLatency: %.2f ms\nJitter: %.2f ms', ...
        firstSpikeReliability, firstSpikeLatency * 1000, firstSpikeJitter * 1000), ...
        'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');

        % Before we extract the waveforms, I want to extract additional
        % waveforms that are not around the optoStim. We can achieve this
        % by looking at the optoDuration and taking some spikes outside of
        % the range.
        firstOptoStim = optoDuration.optoTimes(1);
        lastOptoStim = optoDuration.optoTimes(end);
        temp_spikes = find(unitSpikes < firstOptoStim | unitSpikes > lastOptoStim);
        if length(firstOptoWaveformNumber) < length(temp_spikes)
            spikesToAdd = randsample(temp_spikes,length(firstOptoWaveformNumber));
        else
            spikesToAdd = randsample(temp_spikes,min(120,round(length(temp_spikes)/2)));
        end
        nonOptoWaveformNumber = sort([nonOptoWaveformNumber spikesToAdd']);
        
        % we remove all Nans
        nonOptoWaveformNumber = nonOptoWaveformNumber(~isnan(nonOptoWaveformNumber));
        firstOptoWaveformNumber = firstOptoWaveformNumber(~isnan(firstOptoWaveformNumber));
        otherOptoWaveformNumber = otherOptoWaveformNumber(~isnan(otherOptoWaveformNumber));

        % Now for a given condition and the waveforms that I have
        % extracted, I want to grab the stimulus waveforms.
        myKsDir = cd;
        bindataDir = fullfile(cd, '..');
        gwfparams.bindataDir = bindataDir;
        gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
        apD = dir(fullfile(bindataDir, '*ap*.bin')); % AP band file from spikeGLX specifically
        gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
        gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
        gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
        gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
        gwfparams.nWf = length(firstOptoWaveformNumber);                    % Number of waveforms per unit to pull out
        gwfparams.spikeTimesSamples = sp.ss(sp.clu == goodClusters(i));
        gwfparams.spikeTimesSeconds = sp.st(sp.clu == goodClusters(i));
        gwfparams.spikeClusters = sp.clu(sp.clu==goodClusters(i));
        

        % Now we extract the corresponding waveforms
        wf_first = getThisWaveform(gwfparams,firstOptoWaveformNumber);
        wf_other = getThisWaveform(gwfparams,otherOptoWaveformNumber);
        wf_nonOpto = getThisWaveform(gwfparams,nonOptoWaveformNumber);
    
        % Just make sure to use the maxChannel from wf_nonOpto (because
        % this is the most accurate representation)
        subplot(length(optoConds),4,subplotIdx+1) % Plot the waveforms across trials  
        for wfi = 1:length(nonOptoWaveformNumber)
            plot(1:82,squeeze(wf_nonOpto.waveForms(wfi,wf_nonOpto.maxChannel,:)) - squeeze(wf_nonOpto.waveForms(wfi,wf_nonOpto.maxChannel,1)),...
                '-','Color',[0 0 0 0.2])
            hold on
        end

        for wfi = 1:length(firstOptoWaveformNumber)
            plot(84:165,squeeze(wf_first.waveForms(wfi,wf_nonOpto.maxChannel,:)) - squeeze(wf_first.waveForms(wfi,wf_nonOpto.maxChannel,1)),...
                '-','Color',[0 0 0 0.2])
            hold on
        end

        for wfi = 1:length(otherOptoWaveformNumber)
            plot(250:331,squeeze(wf_other.waveForms(wfi,wf_nonOpto.maxChannel,:)) - squeeze(wf_other.waveForms(wfi,wf_nonOpto.maxChannel,1)),...
                '-','Color',[0 0 0 0.2])
            hold on
        end
        
        xlim([0 333])
        title('Non-opto,first,and other waveforms')

        subplot(length(optoConds),4,subplotIdx+2)
        nonOptoWFMean = wf_nonOpto.waveFormsMean(wf_nonOpto.maxChannel,:) - wf_nonOpto.waveFormsMean(wf_nonOpto.maxChannel,1);
        plot(nonOptoWFMean,'k-','LineWidth',1.5) 
        hold on;
        firstWFMean = wf_first.waveFormsMean(wf_nonOpto.maxChannel,:) - wf_first.waveFormsMean(wf_nonOpto.maxChannel,1);
        plot(firstWFMean,'r-','LineWidth',1.5) 
        hold on;
        otherWFMean = wf_other.waveFormsMean(wf_nonOpto.maxChannel,:) - wf_other.waveFormsMean(wf_nonOpto.maxChannel,1);
        plot(otherWFMean,'b-','LineWidth',1.5) 
    
        r = corrcoef([nonOptoWFMean' firstWFMean' otherWFMean']);
        isupper = logical(triu(ones(size(r)),1));
        r(isupper) = NaN;

        subplot(length(optoConds),4,subplotIdx+3)
        hold off;
        h = heatmap(r,'MissingDataColor','w');
        labels = ["non-top","first wf","other Opto wf"];
        h.XDisplayLabels = labels;
        h.YDisplayLabels = labels; 
        ax = gca;
        ax.FontSize = 10;
    end

    sgtitle(['Cluster ',num2str(goodClusters(i)+1)])
    fPath = fullfile(myKsDir,'OptoScatter');
    try
        saveFigurePdf(unitSpikeOptoScatter, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))    
    catch ME
        disp(['failed to save figure'])
    end

    close(unitSpikeOptoScatter)
end

end



