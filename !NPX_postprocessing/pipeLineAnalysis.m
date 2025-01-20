clear all
clc
format longg 

%spikeTimes = readNPY('spike_times.npy');
%spikeTimes_seconds = readmatrix('spike_times.txt');
myKsDir = cd; % cd to the kilosort output folder
Ksparams.excludeNoise = false; % we make sure all the spike times are imported. 
sp = loadKSdir(myKsDir,Ksparams); % we load the structure in.
% when we load the structure using loadKSdir, spike times is wrong (because
% it assumes 30000Hz sampling rate). Re-import the timestamps (with code
% below)
spiketimes = readmatrix('spike_times.txt'); 
sp.st = spiketimes; 
sp.ss = readNPY('spike_times.npy'); % in samples


% [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
% templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps); 

goodClu_idx = find(sp.cgs == 2); % find index of good clusters
goodClusters = sp.cids(goodClu_idx); % this is -1 from the MATLAB indexing

%% Run quality metric calculation script
[clusterIDs, unitQuality, contaminationRate, LR] = sqKilosort.maskedClusterQualityKilosort(myKsDir);

load("qualMatrix_20.mat")

%% Load the good clusters and plot its waveform
% Here we also grab their:
% 1. ISI violations (in %, less than 2ms?)
% 2. Number of spikes
% 3. L-Ratio
% 4. Waveform from C_waves

%%%%%%%%%%% First calculate the number of ISI violations
isiV = sqKilosort.isiViolations(myKsDir);

%%%%%%%%%%% Then read the waveforms in from C_Waves output
cwavesPath = [cd '\C_waves\'];
c_wf = readNPY([cwavesPath 'mean_waveforms.npy']); % read in the mean waveforms from Cwaves
%snr = readNPY([cwavesPath 'cluster_snr.npy']);
qMetric_vals = [];

for i = 1:length(goodClusters)
    % Load in relevant spike sort quality from isiV
    idx = find(isiV(:,1) == goodClusters(i)+1);
    lenOfSpikeTrain = isiV(idx,4); % this is number of spikes in the cluster
    percISIViolations = isiV(idx,3); % this is percentage of ISI violations (2ms)

    % Below I am defining the parameters for getting the waveforms
    bindataDir = fullfile(cd, '..');
    gwfparams.bindataDir = bindataDir;
    gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
    apD = dir(fullfile(bindataDir, '*ap*.bin')); % AP band file from spikeGLX specifically
    gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
    gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
    gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
    gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf = 1000;                    % Number of waveforms per unit to pull out
    gwfparams.spikeTimesSamples = sp.ss(sp.clu==goodClusters(i)); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeTimesSeconds = sp.st(sp.clu==goodClusters(i));
    gwfparams.spikeClusters = sp.clu(sp.clu==goodClusters(i));
    gwfparams.lenOfSpikeTrain = lenOfSpikeTrain;

    % Rescale the waveforms into uV scale
    wf = getWaveForms_CA(gwfparams);
    squeezedWf = (squeeze(wf.waveFormsMean) * 0.5 / 8192) / 80 * 1000000;     
    squeezedWfStart = (squeeze(wf.waveFormsStartMean) * 0.5 / 8192) / 80 * 1000000;    
    squeezedWfMid = (squeeze(wf.waveFormsMidMean) * 0.5 / 8192) / 80 * 1000000;    
    squeezedWfEnd = (squeeze(wf.waveFormsEndMean) * 0.5 / 8192) / 80 * 1000000;    
    
    % Find the channel for the current unit
    totalWaveform = [];
    zerod_squeezedWf = squeezedWf - squeezedWf(:,1);
    for j = 1:length(sp.xcoords)
        % Iterate through the channels and then find min peak and maximum
        % peak for all of them. Where the maximum peak comes after the
        % minimum
        [minVal minCol] = min(zerod_squeezedWf(j,1:end)); 
        [maxVal maxCol] = max(zerod_squeezedWf(j,minCol:end));
        totalWaveform(j,1) = abs(minVal);       % Find amplitude of the spike
        totalWaveform(j,2) = (maxCol - 1) * 1/30000 * 1000;      % Because maxCol index starts from the index of minCol (this is waveform width)
        % This line is supposed to be (maxCol + minCol) - minCol - 1
        
        peakOfFirst = zerod_squeezedWf(j,1);      % Take the first value of the trace, this should be an okay approx for baseline   
        [firstPeakVal firstPeakCol] = max(zerod_squeezedWf(j,1:minCol));      % Calculate the value of first peak (a)
        [secondPeakVal secondPeakCol] = max(zerod_squeezedWf(j,minCol:end));  % Calculate the value of second peak (b)
        
        a = firstPeakVal - peakOfFirst;     
        b = secondPeakVal - peakOfFirst;           
        totalWaveform(j,3) = (b-a)/(a+b);       % Symmetry index (here known as spikePeakRatio); See Sirota et al., 2008
    end
    [maxAmplitudeWaveform maxChannel] = max(totalWaveform);     % Take the waveform with maximum amplitude and grab the maximum channel 
    spikeWidthTime = totalWaveform(maxChannel(1),2);            % Grab the width of the spike and the ratio from the max channel 
    spikePeakRatio = totalWaveform(maxChannel(1),3);
    rawWf = (squeeze(wf.waveForms) * 0.5 / 8192) / 80 * 1000000;
    rawWf = squeeze(rawWf(:,maxChannel(1),:));

    % Load in relevant C_waves waveforms 
    cwaves_waveformPlot = squeeze(c_wf(goodClusters(i)+1,maxChannel(1),:));

    % To extract spike amplitude, take the maxChannel(1), go through
    % waveforms and then extract
    spikeTimes = sp.st(sp.clu == goodClusters(i)); % already in seconds

    gsaparams.bindataDir = bindataDir;
    gsaparams.dataDir = myKsDir;    % KiloSort/Phy output folder
    apD = dir(fullfile(bindataDir, '*ap*.bin')); % AP band file from spikeGLX specifically
    gsaparams.fileName = apD(1).name;         % .dat file containing the raw 
    gsaparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
    gsaparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
    gsaparams.spikeTimesSamples = sp.ss(sp.clu==goodClusters(i)); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gsaparams.spikeClusters = sp.clu(sp.clu==goodClusters(i));
    gsaparams.lenOfSpikeTrain = lenOfSpikeTrain;
    gsaparams.maxChannel = maxChannel(1);
    gsaparams.nCh = 385;

    sa = getSpikeAmplitudes(gsaparams);

    % Calculate spike rate across time
    binSize = 20; % 10 second size bins
    bins = 0:binSize:max(spikeTimes)+binSize;
    spikesInBins = histcounts(spikeTimes,bins);
    spikesInBins(end) = spikesInBins(end) - 1; % this padding is important to remove the last spike.
    spikesInHz = histcounts(spikeTimes,bins)/binSize; % convert to Hz
    
    spikesAmplitudeInBins = cell(length(spikesInBins), 1);
    current_index = 1;
    for l = 1:length(spikesInBins) % sa.amplitude will have 1 less spike than the amplitude
        split_size = spikesInBins(l);
        spikesAmplitudeInBins{l} = sa.amplitude(current_index:current_index+split_size-1);
        current_index = current_index + split_size;
    end

    %%%%%%%%%%%%%%%% Plot the relevant figures %%%%%%%%%%%%%%%%
    waveFormFigure  = figure('Renderer', 'painters', 'Position', [10 10 1900 1080]);

    % Plot the waveform at the max amplitude channel
    subplot(2,4,1)
    for p = 1:1000
        plot(rawWf(p,:)-rawWf(p,1),'Color',[0 0 0 0.08],'LineWidth',0.5)
        hold on;
    end
    plot(squeezedWf(maxChannel(1),:)-squeezedWf(maxChannel(1),1),'r','LineWidth',2.0) % at peak channel
    set(gca,'box','off')
    set(gca,'FontSize',14)
    xticks([0 15 30 45 60 75 90])         % Points are still in samples
    xticklabels({'0','0.5','1','1.5','2','2.5','3.0'})
    xlabel('Time (ms)')  
    ylabel('\muV')
    title(['Channel ', num2str(maxChannel(1)),...
        ' SpikeWidth = ', num2str(spikeWidthTime),' SpikeAmp = ', num2str(totalWaveform(maxChannel(1),1))])
    minYlim = min(squeezedWf(maxChannel(1),:)-squeezedWf(maxChannel(1),1));
    maxYlim = max(squeezedWf(maxChannel(1),:)-squeezedWf(maxChannel(1),1));
    ylim([minYlim + minYlim/2 maxYlim + 2*maxYlim])

    % Plot waveform across multiple channels
    subplot(2,4,2)
    if maxChannel(1) <= 16
        for nCh = 0:15
            if rem(abs(nCh),2) == 0 % If even channels up or down, plot on the same scale
                plot(1:82,squeezedWf(maxChannel(1)+nCh,:)-squeezedWf(maxChannel(1)+nCh,1)+nCh*3,'k','LineWidth',2.0) % at peak channel
                hold on; 
            else % If odd channels up or down, we should displace the plot
                plot(83:164,squeezedWf(maxChannel(1)+nCh,:)-squeezedWf(maxChannel(1)+nCh,1)+nCh*3+10,'k','LineWidth',2.0) % at peak channel
                hold on; 
            end
        end
    elseif maxChannel(1) >= 367
        for nCh = -15:0
            if rem(abs(nCh),2) == 0 % If even channels up or down, plot on the same scale
                plot(1:82,squeezedWf(maxChannel(1)+nCh,:)-squeezedWf(maxChannel(1)+nCh,1)+nCh*3,'k','LineWidth',2.0) % at peak channel
                hold on; 
            else % If odd channels up or down, we should displace the plot
                plot(83:164,squeezedWf(maxChannel(1)+nCh,:)-squeezedWf(maxChannel(1)+nCh,1)+nCh*3+10,'k','LineWidth',2.0) % at peak channel
                hold on; 
            end
        end
    else
        for nCh = -15:15
            if rem(abs(nCh),2) == 0 % If even channels up or down, plot on the same scale
                plot(1:82,squeezedWf(maxChannel(1)+nCh,:)-squeezedWf(maxChannel(1)+nCh,1)+nCh*3,'k','LineWidth',2.0) % at peak channel
                hold on; 
            else % If odd channels up or down, we should displace the plot
                plot(83:164,squeezedWf(maxChannel(1)+nCh,:)-squeezedWf(maxChannel(1)+nCh,1)+nCh*3+10,'k','LineWidth',2.0) % at peak channel
                hold on; 
            end
        end
    end
    plot(squeezedWf(maxChannel(1),:)-squeezedWf(maxChannel(1),1),'r','LineWidth',2.0) % at peak channel
    set(gca,'box','off')
    set(gca,'XColor', 'none','YColor','none')
    xlim([0 164])

    % Plot the waveform extracted from C_waves and waveforms from
    % beginning, mid, and end
    subplot(2,4,3)
    plot(1:82,cwaves_waveformPlot-cwaves_waveformPlot(1),'k','LineWidth',2.0)
    hold on;
    plot(84:165,squeezedWfStart(maxChannel(1),:)-squeezedWfStart(maxChannel(1),1),'k','LineWidth',2.0)
    hold on;
    plot(167:248,squeezedWfMid(maxChannel(1),:)-squeezedWfMid(maxChannel(1),1),'r','LineWidth',2.0)
    hold on;
    plot(250:331,squeezedWfEnd(maxChannel(1),:)-squeezedWfEnd(maxChannel(1),1),'b','LineWidth',2.0)
    hold on;
    set(gca,'box','off')
    set(gca,'FontSize',14)
    set(gca,'XColor', 'none','YColor','none')
    ylabel('\muV')
    minYlim = min(squeezedWf(maxChannel(1),:)-squeezedWf(maxChannel(1),1));
    maxYlim = max(squeezedWf(maxChannel(1),:)-squeezedWf(maxChannel(1),1));
    ylim([minYlim + minYlim/2 maxYlim + 2*maxYlim])
    lgd1 = legend('Cwaves','Beginning','Mid','End','Location','Southeastoutside');
    lgd1.FontSize = 9;

    % Plot spike amplitude and spike rate over time
    subplot(4,3,[7 8])
    histogram('BinEdges',bins,'BinCounts',spikesInHz,'FaceColor','black','EdgeColor','none') % spike rate over time
    ylabel('Spike rate (Hz)')
    set(gca,'XColor', 'none')
    xlim([0 max(spikeTimes)])
    set(gca,'FontSize',10)
    set(gca,'box','off')

    % Add amplitude across time
    subplot(4,3,[10 11])
    nonzeroAmp = find(sa.amplitude/mean(sa.amplitude) ~= 0);
    scatter(spikeTimes(nonzeroAmp),nonzeros(sa.amplitude/mean(sa.amplitude)),'filled')
    ylabel('Amplitude normalized to mean')
    ylim([0 3])
    xlabel('Time (s)')
    xlim([0 max(spikeTimes)])
    set(gca,'FontSize',10)
    set(gca,'box','off')

    % Plot amplitude histogram (but normalized to the mean here)
    subplot(4,6,[23])
    nonzeroAmp = find(sa.amplitude/mean(sa.amplitude) ~= 0);
    histogram(nonzeros(sa.amplitude/mean(sa.amplitude)),'orientation','horizontal','Normalization','probability')
    xlabel('Proportion')
    ylim([0 3])
    set(gca,'FontSize',10)
    set(gca,'box','off')

    % Here I am going to generate a separate plot that plots waveforms and
    % its correlation across every 5 minutes
    waveformToPlot = nan(size(wf.waveFormsMinMean,1),82);
    for p = 1:size(wf.waveFormsMinMean,1)
        waveformToPlot(p,:) = squeeze(wf.waveFormsMinMean(p,maxChannel(1),:)) - squeeze(wf.waveFormsMinMean(p,maxChannel(1),1));  
    end 
    waveformCrossCoef=corrcoef(waveformToPlot'); % This calculates cross coef between two waveforms
    isupper = logical(triu(ones(size(waveformCrossCoef)),1));
    waveformCrossCoef(isupper) = NaN;
    checkWaveformCrossCoef = diag(waveformCrossCoef,-1);
    
    subplot(4,6,17)
    corrcoefBool = 0;
    h = heatmap(waveformCrossCoef,'MissingDataColor','w','ColorLimits',[0 1],'Colormap',parula);
    if isempty(find(checkWaveformCrossCoef <= 0.9)) & isempty(find(waveformCrossCoef <= 0.7))
        title('Good corrcoef')
        corrcoefBool = 1;
    else
        title('Bad corrcoef')
        corrcoefBool = 0;
    end
    
    subplot(2,6,12)
    histogram(diff(spikeTimes)*1000,0:1:100,'FaceColor','black','Edgecolor','none',...
        'Normalization','probability')
    hold on;
    xline(2,'Color','r','LineWidth',2)
    set(gca,'box','off')
    set(gca,'FontSize',10)
    xlabel('ISI (ms)')
    xlim([0 20])

    subplot(2,4,4)
    scatter(sp.xcoords,sp.ycoords,'ks','filled')
    hold on
    scatter(sp.xcoords(maxChannel(1)),sp.ycoords(maxChannel(1)),200,'rs','filled')
    xlabel('\mum')
    ylabel('From deepest part of probe')

    LRatioUnit = qualMatrix(find(qualMatrix(:,1) == goodClusters(i)+1),4);

    sgtitle(['Cluster ',num2str(goodClusters(i)+1),'; NumOfSpikes= ',num2str(lenOfSpikeTrain),...
        '; % ISI violations= ',num2str(percISIViolations),'% ; L-Ratio= ',num2str(LRatioUnit)])

    % Save the figures
    fPath = fullfile(myKsDir,'waveforms');
    try
        saveFigurePdf(waveFormFigure, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
    catch ME
        disp(['failed to save figure'])
    end
    close(waveFormFigure);

    % Save the variables
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).amplitude = totalWaveform(maxChannel(1),1); 
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).maxChannel = maxChannel(1);
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).spikeWidth = spikeWidthTime;
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).numOfSpikes = lenOfSpikeTrain;
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).ISIv = percISIViolations;
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).LRatio = LRatioUnit;
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).waveformCorrelation = corrcoefBool;
    qMetric_vals.(sprintf('Unit%d',goodClusters(i)+1)).cluster = goodClusters(i)+1;

end
close all

% Save qMetric_vals structure
fPath = fullfile(myKsDir,'qMetrics');

% Ensure the directory exists (optional)
if ~isfolder(fPath)
    mkdir(fPath); % Create the directory if it doesn't exist
end

save(fullfile(fPath,'qMetric_vals.mat'),'qMetric_vals')

