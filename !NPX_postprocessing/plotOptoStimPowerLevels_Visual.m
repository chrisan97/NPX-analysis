function plotOptoStimPowerLevels_Visual(sp, goodClusters, visStimTimes, myKsDir, mw_Table)
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 9/2/2024
%
% ------
% Inputs
% ------
% sp: spike structure holding array of information. Usually generated from
% loadKsDir
% goodClusters: array holding cluster numbers that are labeled 'good'
% visStimTimes: visual stimulus timing in array
% myKsDir: directory to the kilosort folder
% mw_Table: usually n x 2 table where 1st column is the power and the
% second column is the correpsonding time of optogenetic stimulation
%
% ------
% Outputs
% ------
% Generates a plot for each unit, use this when you have optogenetics
%


%% First, extract the power levels in this recording
powerLevels = unique(mw_Table.Power);
binSize = 0.050;
binSize_zoom = 0.005;
edges = -0.7:binSize:3.0;
edgesZoom = -0.1:binSize_zoom:0.1;
numOfbins = length(edges)-1; % calculate the number of bins
numOfBins_zoom = length(edgesZoom)-1;

colors = [0 0 0; 1 0 1; 0 1 1; 1 0.7 0; 0.3 0 0.7]; % black, magenta, cyan, etc


%% We now want to iterate through each power levels
for i = 1:length(goodClusters) % for all the good clusters in my recording
    unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds, take spikes from the good cluster
    numOfSpikesPSTH = zeros(length(powerLevels),numOfbins); % initialize the array to hold PSTHs
    numOfSpikesPSTH_zoom = zeros(length(powerLevels),numOfBins_zoom); % initialize the array to hold PSTHs

    unitSpikePowerScatter = figure('Renderer', 'painters', 'Position', [10 10 1200 1200]);
    subplot(2,2,1)
    % Plot across power levels
    y = 0; % for plotting
    numOfStims = height(mw_Table);
    for plvl = 1:length(powerLevels)
        % find matching idx and extract corresponding optoTimes
        thisPowerIdx = find(mw_Table.Power == powerLevels(plvl));
        thisPowerStimTimes = mw_Table.OptoTimes(thisPowerIdx);
        % We now iterate through the stimulation times
        firstOptoSpikeLatency = [];
        for stim = 1:length(thisPowerStimTimes)
            alignSpikes = unitSpikes - thisPowerStimTimes(stim); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, 'k', 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3);
            N = histcounts(indivTrial, edges);
            numOfSpikesPSTH(plvl,:) = numOfSpikesPSTH(plvl,:) + N;

            % Calculate latency to first spike
            firstOptoSpikeLatency = [firstOptoSpikeLatency indivTrial(find(indivTrial >= 0,1,'first'))];
        end
       
    % draw a red line at end of every plvl
    yline(y,'r-','LineWidth',1.5)
    end
    % Plot other aspects of the scatter (line at 0, etc)
    plot(zeros(1, numOfStims), 1:numOfStims, 'k--', 'LineWidth', 1.0)
    hold on;
    fill([0, 0.5, 0.5, 0], [0, 0, numOfStims, numOfStims], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
    xlim([-0.7, 3.0]);
    ylim([0, numOfStims]);
    xlabel('Time (s)');
    ylabel('Optogenetics pulse');
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.FontSize = 14;
    box off;
    set(gca,'LineWidth',2,'FontSize',14)

    subplot(2,2,2)
    % Plot across power levels
    y = 0; % for plotting
    numOfStims = height(mw_Table);
    for plvl = 1:length(powerLevels)
        % find matching idx and extract corresponding optoTimes
        thisPowerIdx = find(mw_Table.Power == powerLevels(plvl));
        thisPowerStimTimes = mw_Table.OptoTimes(thisPowerIdx);
        % We now iterate through the stimulation times
        firstOptoSpikeLatency = [];
        for stim = 1:length(thisPowerStimTimes)
            alignSpikes = unitSpikes - thisPowerStimTimes(stim); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, 'k', 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            indivTrial = alignSpikes(alignSpikes < 0.1 & alignSpikes >= -0.1);
            N = histcounts(indivTrial, edgesZoom);
            numOfSpikesPSTH_zoom(plvl,:) = numOfSpikesPSTH_zoom(plvl,:) + N;
        end
    % draw a red line at end of every plvl
    yline(y,'r-','LineWidth',1.5)
    end
    % Plot other aspects of the scatter (line at 0, etc)
    plot(zeros(1, numOfStims), 1:numOfStims, 'k--', 'LineWidth', 1.0)
    hold on;
    fill([0, 0.5, 0.5, 0], [0, 0, numOfStims, numOfStims], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
    xlim([-0.1, 0.1]);
    ylim([0, numOfStims]);
    xlabel('Time (s)');
    ylabel('Optogenetics pulse');
    ax = gca;
    ax.XAxis.Exponent = 0;
    box off;
    set(gca,'LineWidth',2,'FontSize',14)

    %% Plot PSTHs
    subplot(2,2,3)
    for plvl = 1:length(powerLevels)
        thisPowerIdx = find(mw_Table.Power == powerLevels(plvl));
        spikeRatePSTH = numOfSpikesPSTH(plvl,:)/(length(thisPowerIdx) * binSize);
        plot(spikeRatePSTH,'Color',[0 0 1 plvl*0.17],'LineWidth',1.7)
        hold on;
    end
    timePoints = [-0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
    ind = round((timePoints + 0.7)/binSize) + 1;
    xticks(ind);
    xticklabels(timePoints);
    xline(ind(2),'r-','LineWidth',1.2)
    xlabel('Time (s)')
    ylabel('Hz')
    legend(num2str(powerLevels))
    ax = gca;
    box off;
    set(gca,'LineWidth',2,'FontSize',14)

    %% Plot PSTH for zoom
    subplot(2,2,4)
    for plvl = 1:length(powerLevels)
        thisPowerIdx = find(mw_Table.Power == powerLevels(plvl));
        spikeRatePSTH = numOfSpikesPSTH_zoom(plvl,:)/(length(thisPowerIdx) * binSize_zoom);
        plot(spikeRatePSTH,'Color',[0 0 1 plvl*0.17],'LineWidth',1.7)
        hold on;
    end
    timePoints = [-0.1, 0, 0.02, 0.04, 0.06 0.08, 0.1];
    ind = round((timePoints + 0.1)/binSize_zoom) + 1;
    xticks(ind);
    xticklabels(timePoints);
    xline(ind(2),'r-','LineWidth',1.2)
    xlabel('Time (s)')
    ylabel('Hz')
    %legend(num2str(powerLevels))
    ax = gca;
    box off;
    set(gca,'LineWidth',2,'FontSize',14)
    sgtitle(['Cluster ',num2str(goodClusters(i)+1)])

    fPath = fullfile(myKsDir,'PowerLevels');
    try
        saveFigurePdf(unitSpikePowerScatter, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
    catch ME
        disp(['failed to save figure'])
    end

    close(unitSpikePowerScatter)
end



end