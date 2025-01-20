 function plotVisStimOptoScatter(sp, goodClusters, visStimTimes, myKsDir, RPiVisAllMap)
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
%
% ------
% Outputs
% ------
% Generates a plot for each unit, use this when you have optogenetics
%

%%%% Preset colours %%%%
colors = {'k','k','k'};
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Preset binSize %%%%
binSize = 0.020;
edges = -0.3:binSize:2.0;
edgesOrientation = 0:binSize:1.5;
edgesPartial = 0.5:binSize:1.0;
%%%%%%%%%%%%%%%%%%%%%%%%

assert(mod(length(visStimTimes),120) == 0) % Check that visual stimulation time is a multiple of 120

numOfbins = length(edges)-1; % calculate the number of bins
visMapConds = unique(RPiVisAllMap.Var1); % These are list of visual mappings that the animal viewed in this session
visMapIdentity = unique(RPiVisAllMap.Var2); % This should be an array holding visual stimulation angesl
optoIdentity = []; % this should hold the type of optogenetics stimulation that was delivered alongside the visual stimulus
% optoIdentity = -1 means no optogenetics interleaved
% optoIdentity = 0 means no optogenetics but interleaved
% optoIdentity = 1 means yes optogenetics and not interleaved
for i = 1:length(RPiVisAllMap.Var3)
    temp = RPiVisAllMap.Var3(i);
    if isnan(temp)
        optoIdentity = [optoIdentity -1];
        RPiVisAllMap.Var3(i) = -1;
    else
        optoIdentity = [optoIdentity temp];
    end
end

for i = 1:length(goodClusters)
    unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds, take spikes from the good cluster
    numOfSpikesPSTH = zeros(length(visMapConds),numOfbins); % initialize the array to hold PSTHs
    numOfSpikesPSTH_noOpto = zeros(length(visMapConds),numOfbins); % no opto
    numOfSpikesPSTH_opto = zeros(length(visMapConds),numOfbins); % to hold opto

    numOfSpikesOrientation = {};
    numOfSpikesOrientation_partial = {};
    numOfSpikesOrientation_noOpto = {};
    numOfSpikesOrientation_partial_noOpto = {};
    numOfSpikesOrientation_opto = {};
    numOfSpikesOrientation_partial_opto = {};

    unitSpikePSTHScatter = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);

    %% First we want to plot the scatter for the grating conditions
    for cond = 1:length(visMapConds)
        rowIdx = mod(cond-1,length(visMapConds)) + 1;
        colIdx = ceil(cond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;

        % Initialize subplot to plot scatters. 
        subplot(length(visMapConds),4,subplotIdx)
        y = 0; % for plotting
        numOfStims = sum(matches(RPiVisAllMap.Var1,visMapConds(cond))); % how many instances of the stimulation in this condition
        thisVisStimTimes = RPiVisAllMap.StimTimes(matches(RPiVisAllMap.Var1,visMapConds(cond)));
        assert(length(thisVisStimTimes) == numOfStims) % just make sure extraction was correct

        for stim = 1:numOfStims % for each stimulation
            alignSpikes = unitSpikes - thisVisStimTimes(stim); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3);
            N = histcounts(indivTrial, edges);
            if optoIdentity(stim) == 0
                numOfSpikesPSTH_noOpto(cond,:) = numOfSpikesPSTH_noOpto(cond,:) + N; % no stimulation interleaved
            elseif optoIdentity(stim) == 1
                numOfSpikesPSTH_opto(cond,:) = numOfSpikesPSTH_opto(cond,:) + N; % stimulation interleaved
            else 
                numOfSpikesPSTH(cond,:) = numOfSpikesPSTH(cond,:) + N; % no stimulation block
            end
        end
        % Plot other aspects of the scatter (line at 0, etc)
        plot(zeros(1, numOfStims+1), 1:numOfStims+1, 'k--', 'LineWidth', 1.0)
        hold on;
        fill([0, 1.5, 1.5, 0], [0, 0, numOfStims, numOfStims], 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none')
        xlim([-0.3, 2.0]);
        ylim([0, numOfStims+1]);
        xlabel('Time (s)');
        ylabel('Ordered by trial order');
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10;
    end
    
    %% Second, now plotting the corresponding PSTHs
    for cond = 1:length(visMapConds)
        plotCond = cond + length(visMapConds);
        rowIdx = mod(plotCond-1,length(visMapConds)) + 1;
        colIdx = ceil(plotCond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;

        subplot(length(visMapConds),4,subplotIdx)
        numOfStims = sum(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == -1)); % how many instances of the stimulation in this condition
        numOfStims_noOpto = sum(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 0));
        numOfStims_opto = sum(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 1));

        spikeRatePSTH = numOfSpikesPSTH(cond,:) / (numOfStims * binSize);
        spikeRatePSTH_noOpto = numOfSpikesPSTH_noOpto(cond,:) / (numOfStims_noOpto * binSize);
        spikeRatePSTH_opto = numOfSpikesPSTH_opto(cond,:) / (numOfStims_opto * binSize);

        combinedTemp = cat(1,spikeRatePSTH,spikeRatePSTH_noOpto,spikeRatePSTH_opto);
        maxYlim = max(combinedTemp,[],'all');

        plot(spikeRatePSTH,'Color',[0 0 0],'LineWidth',1.2)
        hold on;
        plot(spikeRatePSTH_noOpto,'Color',[1 0 0],'LineWidth',1.2)
        hold on;
        plot(spikeRatePSTH_opto,'Color',[0 0 1],'LineWidth',1.2)

        % Plot other aspects of the PSTH 
        timePoints = [-0.3, 0, 0.5, 1.0, 1.5, 2.0];
        ind = round((timePoints + 0.3)/binSize) + 1;
        xticks(ind);
        xticklabels(timePoints);
        hold on;
        fill([ind(2), ind(5), ind(5), ind(2)], [0, 0, maxYlim + 3, maxYlim + 3], 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none')
        hold on;
        fill([ind(3), ind(4), ind(4), ind(3)], [0, 0, maxYlim + 3, maxYlim + 3], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
        ylim([0, maxYlim + 3]);
        xlabel('Time (s)');
        ylabel('Hz');
        legend('Pure orientation','No Opto','Opto')
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10;
    end

    %% Then, we plot the same scatters with PSTH but organized by the order of the visual stimulus presentation
    for cond = 1:length(visMapConds)
        plotCond = cond + 2*length(visMapConds);
        rowIdx = mod(plotCond-1,length(visMapConds)) + 1;
        colIdx = ceil(plotCond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;
        
        % Initialize the subplot for ordered scattered
        subplot(length(visMapConds),4,subplotIdx)
        numOfStims = sum(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == -1 )); % how many instances of the stimulation in this condition
        numOfStims_noOpto = sum(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 0));
        numOfStims_opto = sum(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 1));

        % We are now going to get visual stimulus times for each condition
        thisVisStimTimes = RPiVisAllMap.StimTimes(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == -1));
        thisVisStimOrder = RPiVisAllMap.Var2(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == -1));
        thisVisStimTimes_noOpto = RPiVisAllMap.StimTimes(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 0));
        thisVisStimOrder_noOpto = RPiVisAllMap.Var2(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 0));
        thisVisStimTimes_opto = RPiVisAllMap.StimTimes(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 1));
        thisVisStimOrder_opto = RPiVisAllMap.Var2(matches(RPiVisAllMap.Var1,visMapConds(cond)) & (RPiVisAllMap.Var3 == 1));
        
        % This is a N by 2 matrix, where first row holds the identity
        % of the grating and the second row holds the timestamps
        sortedThisVisStimOrder = sortrows([thisVisStimOrder thisVisStimTimes]); 
        sortedThisVisStimOrder_noOpto = sortrows([thisVisStimOrder_noOpto thisVisStimTimes_noOpto]); 
        sortedThisVisStimOrder_opto = sortrows([thisVisStimOrder_opto thisVisStimTimes_opto]); 

        assert(length(sortedThisVisStimOrder) == numOfStims) % just make sure extraction was correct
        assert(length(sortedThisVisStimOrder_noOpto) == numOfStims_noOpto) % just make sure extraction was correct
        assert(length(sortedThisVisStimOrder_opto) == numOfStims_opto) % just make sure extraction was correct

        temp_fullDuration = [];
        temp_partialDuration = [];
        temp_fullDuration_noOpto = [];
        temp_partialDuration_noOpto = [];
        temp_fullDuration_opto = [];
        temp_partialDuration_opto = [];
        y = 0;
        for stim = 1:length(sortedThisVisStimOrder) % for each stimulation
            alignSpikes = unitSpikes - sortedThisVisStimOrder(stim,2); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            thisVisStimIdentity = sortedThisVisStimOrder(stim,1);
            
            % for full duration
            indivTrial = alignSpikes(alignSpikes <= 1.5 & alignSpikes >= 0);
            N = histcounts(indivTrial, edgesOrientation);
            temp_fullDuration = [temp_fullDuration; N];
            
            % for partial
            indivTrial = alignSpikes(alignSpikes <= 1.0 & alignSpikes >= 0.5);
            N = histcounts(indivTrial, edgesPartial);
            temp_partialDuration = [temp_partialDuration; N];
        end 
        % Draw a horizontal line
        yline(y,'r-','LineWidth',2)
        % For no opto
        for stim = 1:length(sortedThisVisStimOrder_noOpto)
            alignSpikes = unitSpikes - sortedThisVisStimOrder_noOpto(stim,2); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            % for full duration
            indivTrial = alignSpikes(alignSpikes <= 1.5 & alignSpikes >= 0);
            N = histcounts(indivTrial, edgesOrientation);
            temp_fullDuration_noOpto = [temp_fullDuration_noOpto; N];

            % for partial
            indivTrial = alignSpikes(alignSpikes <= 1.0 & alignSpikes >= 0.5);
            N = histcounts(indivTrial, edgesPartial);
            temp_partialDuration_noOpto = [temp_partialDuration_noOpto; N];
        end
        % Draw a horizontal line
        yline(y,'r-','LineWidth',2)
        % For opto trials
        for stim = 1:length(sortedThisVisStimOrder_opto)
            alignSpikes = unitSpikes - sortedThisVisStimOrder_opto(stim,2); % align the spikes to the stimulation onset                scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            y = y + 1;               
            % for full duration
            indivTrial = alignSpikes(alignSpikes <= 1.5 & alignSpikes >= 0);
            N = histcounts(indivTrial, edgesOrientation);
            temp_fullDuration_opto = [temp_fullDuration_opto; N];

            % for partial
            indivTrial = alignSpikes(alignSpikes <= 1.0 & alignSpikes >= 0.5);
            N = histcounts(indivTrial, edgesPartial);
            temp_partialDuration_opto = [temp_partialDuration_opto; N];
        end
        % Plot other aspects of the scatter (line at 0, etc)
        plot(zeros(1, numOfStims+numOfStims_noOpto+numOfStims_opto+1), 1:numOfStims+numOfStims_noOpto+numOfStims_opto+1, 'k--', 'LineWidth', 1.0)
        hold on;
        fill([0, 1.5, 1.5, 0], [0, 0, numOfStims+numOfStims_noOpto+numOfStims_opto+1, numOfStims+numOfStims_noOpto+numOfStims_opto+1], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
        xlim([-0.3, 2.0]);
        ylim([0, numOfStims+numOfStims_noOpto+numOfStims_opto+1]);
        xlabel('Time (s)');
        ylabel('Ordered by orientation');
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10;    

        numOfSpikesOrientation{cond,1} = temp_fullDuration;
        numOfSpikesOrientation_partial{cond,1} = temp_partialDuration;
        numOfSpikesOrientation{cond,2} = sortedThisVisStimOrder(:,1);
        numOfSpikesOrientation_partial{cond,2} = sortedThisVisStimOrder(:,1);

        numOfSpikesOrientation_noOpto{cond,1} = temp_fullDuration_noOpto;
        numOfSpikesOrientation_partial_noOpto{cond,1} = temp_partialDuration_noOpto;
        numOfSpikesOrientation_noOpto{cond,2} = sortedThisVisStimOrder_noOpto(:,1);
        numOfSpikesOrientation_partial_noOpto{cond,2} = sortedThisVisStimOrder_noOpto(:,1);

        numOfSpikesOrientation_opto{cond,1} = temp_fullDuration_opto;
        numOfSpikesOrientation_partial_opto{cond,1} = temp_partialDuration_opto;
        numOfSpikesOrientation_opto{cond,2} = sortedThisVisStimOrder_opto(:,1);
        numOfSpikesOrientation_partial_opto{cond,2} = sortedThisVisStimOrder_opto(:,1);
    end

    
    %% Last, plotting the orientation plots 
    % numOfSpikesOrientation is a cell array, where the first column
    % holds the PSTH for each stimulus presentation
    % second column holds the stimulus order matrix

    sumSpikesOrientation = cellfun(@(x) sum(x,2), numOfSpikesOrientation, 'UniformOutput', false);
    sumSpikesOrientation_partial = cellfun(@(x) sum(x,2), numOfSpikesOrientation_partial, 'UniformOutput', false);
    sumSpikesOrientation_noOpto = cellfun(@(x) sum(x,2), numOfSpikesOrientation_noOpto, 'UniformOutput', false);
    sumSpikesOrientation_partial_noOpto = cellfun(@(x) sum(x,2), numOfSpikesOrientation_partial_noOpto, 'UniformOutput', false);
    sumSpikesOrientation_opto = cellfun(@(x) sum(x,2), numOfSpikesOrientation_opto, 'UniformOutput', false);
    sumSpikesOrientation_partial_opto = cellfun(@(x) sum(x,2), numOfSpikesOrientation_partial_opto, 'UniformOutput', false);

    statsOrientation = {}; % this is a cell array where each row is a condition and column is a visual map identity. 
    statsOrientation_partial = {};
    statsOrientation_noOpto = {}; % this is a cell array where each row is a condition and column is a visual map identity. 
    statsOrientation_partial_noOpto = {};
    statsOrientation_opto = {}; % this is a cell array where each row is a condition and column is a visual map identity. 
    statsOrientation_partial_opto = {};
    % First element of the matrix is the mean and the second element is
    % std, 3rd element is SEM
   
    for cond = 1:length(visMapConds)
        for j = 1:length(visMapIdentity)
            % For full
            temp_ind = find(sumSpikesOrientation{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation{cond,1}(temp_ind)/(max(edgesOrientation)-min(edgesOrientation));
            statsOrientation{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];
        
            temp_ind = find(sumSpikesOrientation_partial{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation_partial{cond,1}(temp_ind)/(max(edgesPartial)-min(edgesPartial));
            statsOrientation_partial{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];     
        
            % For no opto trials
            temp_ind = find(sumSpikesOrientation_noOpto{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation_noOpto{cond,1}(temp_ind)/(max(edgesOrientation)-min(edgesOrientation));
            statsOrientation_noOpto{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];
            
            temp_ind = find(sumSpikesOrientation_partial_noOpto{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation_partial_noOpto{cond,1}(temp_ind)/(max(edgesPartial)-min(edgesPartial));
            statsOrientation_partial_noOpto{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];
            
            % For opto trials
            temp_ind = find(sumSpikesOrientation_opto{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation_opto{cond,1}(temp_ind)/(max(edgesOrientation)-min(edgesOrientation));
            statsOrientation_opto{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];
            
            temp_ind = find(sumSpikesOrientation_partial_opto{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation_partial_opto{cond,1}(temp_ind)/(max(edgesPartial)-min(edgesPartial));
            statsOrientation_partial_opto{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];
        end
    end
    
    % Now plot the orientation stats
    for cond = 1:length(visMapConds)
        plotCond = cond + 3*length(visMapConds);
        rowIdx = mod(plotCond-1,length(visMapConds)) + 1;
        colIdx = ceil(plotCond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;

        subplot(length(visMapConds),4,subplotIdx)
            
        assert(length(statsOrientation) == length(visMapIdentity)) % just make sure extraction was correct
        assert(length(statsOrientation_partial) == length(visMapIdentity)) % just make sure extraction was correct
        assert(length(statsOrientation_noOpto) == length(visMapIdentity)) % just make sure extraction was correct
        assert(length(statsOrientation_partial_noOpto) == length(visMapIdentity)) % just make sure extraction was correct
        assert(length(statsOrientation_opto) == length(visMapIdentity)) % just make sure extraction was correct
        assert(length(statsOrientation_partial_opto) == length(visMapIdentity)) % just make sure extraction was correct

        % all conditions
        fullSEM = cellfun(@(x) x(3), statsOrientation);
        fullMean = cellfun(@(x) x(1), statsOrientation);
        partialSEM = cellfun(@(x) x(3), statsOrientation_partial);
        partialMean = cellfun(@(x) x(1), statsOrientation_partial);

        % no opto 
        fullSEM_noOpto = cellfun(@(x) x(3), statsOrientation_noOpto);
        fullMean_noOpto = cellfun(@(x) x(1), statsOrientation_noOpto);
        partialSEM_noOpto = cellfun(@(x) x(3), statsOrientation_partial_noOpto);
        partialMean_noOpto = cellfun(@(x) x(1), statsOrientation_partial_noOpto);

        % opto only
        fullSEM_opto = cellfun(@(x) x(3), statsOrientation_opto);
        fullMean_opto = cellfun(@(x) x(1), statsOrientation_opto);
        partialSEM_opto = cellfun(@(x) x(3), statsOrientation_partial_opto);
        partialMean_opto = cellfun(@(x) x(1), statsOrientation_partial_opto);

        % errorbar(1:30:360, fullMean(cond,:), fullSEM(cond,:),'k-');
        % hold on;
        % errorbar(1:30:350, fullMean(cond,:), fullSEM(cond,:),'k-','LineWidth',1.4,'Capsize',0);
        % hold on;
        % errorbar(1:30:360, fullMean_noOpto(cond,:), fullSEM_noOpto(cond,:),'r-','LineWidth',1.4,'Capsize',0);
        % hold on;
        % errorbar(1:30:360, fullMean_opto(cond,:), fullSEM_opto(cond,:),'b-','LineWidth',1.4,'Capsize',0);
        % hold on;
        % yyaxis right
        errorbar(1:30:360, partialMean(cond,:), partialSEM(cond,:),'k-','LineWidth',1.4,'Capsize',0);
        hold on;
        errorbar(1:30:360, partialMean_noOpto(cond,:), partialSEM_noOpto(cond,:),'r-','LineWidth',1.4,'Capsize',0);
        hold on;
        errorbar(1:30:360, partialMean_opto(cond,:), partialSEM_opto(cond,:),'b-','LineWidth',1.4,'Capsize',0);

        temp_combined = [partialMean(cond,:) partialMean_noOpto(cond,:) partialMean_opto(cond,:)];

        ylabel('Hz')
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10; 
        xlim([0 360])
        ylim([0 max(temp_combined)+max(temp_combined) * 0.2])
        legend('Pure orientation','No-Opto','Opto','Location','northeast')
    end

    sgtitle(['Cluster ',num2str(goodClusters(i)+1)])
    fPath = fullfile(myKsDir,'OptoOrientation');
    try
        saveFigurePdf(unitSpikePSTHScatter, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
    catch ME
        disp(['failed to save figure'])
    end

    close(unitSpikePSTHScatter);
    orientationStats.cluster =  goodClusters(i);
    orientationStats.fullMean.noOpto = fullMean_noOpto;
    orientationStats.fullMean.opto = fullMean_opto;
    orientationStats.fullMean.pure = fullMean;
    orientationStats.fullSEM.noOpto = fullSEM_noOpto;
    orientationStats.fullSEM.opto = fullSEM_opto;
    orientationStats.fullSEM.pure = fullSEM;
    orientationStats.numSpikes.noOpto = numOfSpikesOrientation_noOpto;
    orientationStats.numSpikes.opto = numOfSpikesOrientation_opto;
    orientationStats.numSpikes.pure = numOfSpikesOrientation;      
    orientationStats.binEdges = edgesOrientation;
    orientationStats.visDirections = visMapIdentity;
    
    fname = append('Cluster_',num2str(goodClusters(i)+1));
    fPath = fullfile(myKsDir,'OptoOrientation','struct',fname);
    save(fPath,'orientationStats')
end
close all;
end