function plotVisStimScatter(sp, goodClusters, myKsDir, RPiVisAllMap)
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/18/2024
% 
% This code is usually used if there are no optogenetics delivered in the
% session. For code where optogenetics has been delivered, check:
% plotVisStimOptoScatter
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
% Generates a plot for each unit
%


%%%% Preset colours %%%%
colors = {'k','k','k'};
%%%%%%%%%%%%%%%%%%%%%%%%
binSize = 0.050;
edges = -0.3:binSize:2.0;
edgesOrientation = 0:binSize:1.5;
edgesPartial = 0.5:binSize:1.5;

numOfbins = length(edges)-1; % calculate the number of bins
visMapConds = unique(RPiVisAllMap.Var1); % These are list of visual mappings that the animal viewed in this session
visMapIdentity = unique(RPiVisAllMap.Var2);

% We iterate the clusters
for i = 1:length(goodClusters)
    unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds, take spikes from the good cluster
    numOfSpikesPSTH = zeros(length(visMapConds),numOfbins); % initialize the array to hold PSTHs
    numOfSpikesOrientation = {};
    numOfSpikesOrientation_partial = {};
    unitSpikePSTHScatter = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);

    %% First we want to plot the scatter for the grating conditions
    for cond = 1:length(visMapConds)
        % We first calculate subplot size and assign the subplot number
        rowIdx = mod(cond-1,length(visMapConds)) + 1;
        colIdx = ceil(cond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;

        % Generate the first subplot for scatters
        subplot(length(visMapConds),4,subplotIdx)
        y = 0; % for plotting
        numOfStims = sum(matches(RPiVisAllMap.Var1,visMapConds(cond))); % how many instances of the stimulation in this condition
        thisVisStimTimes = RPiVisAllMap.StimTimes(matches(RPiVisAllMap.Var1,visMapConds(cond)));
        assert(length(thisVisStimTimes) == numOfStims) % just make sure extraction was correct
        % for each stimulation
        for stim = 1:numOfStims
            alignSpikes = unitSpikes - thisVisStimTimes(stim); % align the spikes to the stimulation onset
            scatter(alignSpikes, ones(1, length(alignSpikes)) + y, 10, colors{cond}, 'filled'); % Plot the scatters
            y = y + 1;
            hold on;
            
            indivTrial = alignSpikes(alignSpikes < max(edges) & alignSpikes >= min(edges)); % Save the spike counts for PSTH
            N = histcounts(indivTrial, edges);
            numOfSpikesPSTH(cond,:) = numOfSpikesPSTH(cond,:) + N;
        end
        % Plot other aspects of the scatter (line at 0, etc)
        plot(zeros(1, numOfStims+1), 1:numOfStims+1, 'k--', 'LineWidth', 1.0)
        hold on;
        fill([0, 1.5, 1.5, 0], [0, 0, numOfStims, numOfStims], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
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
        % We calculate subplot size and assign the subplot number
        plotCond = cond + length(visMapConds);
        rowIdx = mod(plotCond-1,length(visMapConds)) + 1;
        colIdx = ceil(plotCond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;

        % Generate the subplot for PSTH
        subplot(length(visMapConds),4,subplotIdx)
        numOfStims = sum(matches(RPiVisAllMap.Var1,visMapConds(cond))); % how many instances of the stimulation in this condition
        spikeRatePSTH = numOfSpikesPSTH(cond,:) / (numOfStims * binSize); % change spike counts into Hz
        maxYlim = max(spikeRatePSTH);
        plot(spikeRatePSTH)
        % Plot other aspects of the PSTH 
        timePoints = [-0.3, 0, 0.5, 1.0, 1.5, 2.0];   
        ind = round((timePoints + 0.3)/binSize) + 1;
        xticks(ind);
        xticklabels(timePoints);
        hold on;
        fill([ind(2), ind(5), ind(5), ind(2)], [0, 0, maxYlim + 1, maxYlim + 1], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
        ylim([0, maxYlim + 1]);
        xlabel('Time (s)');
        ylabel('Hz');
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10;
    end

    %% Then, we plot the same scatters with PSTH but organized by the order of the visual stimulus presentation
    for cond = 1:length(visMapConds)
        % We calculate subplot size and assign the subplot number
        plotCond = cond + 2*length(visMapConds);
        rowIdx = mod(plotCond-1,length(visMapConds)) + 1;
        colIdx = ceil(plotCond/length(visMapConds));
        subplotIdx = (rowIdx - 1) * 4 + colIdx;
        
        subplot(length(visMapConds),4,subplotIdx)
        y = 0; % for plotting
        numOfStims = sum(matches(RPiVisAllMap.Var1,visMapConds(cond))); % how many instances of the stimulation in this condition
        thisVisStimTimes = RPiVisAllMap.StimTimes(matches(RPiVisAllMap.Var1,visMapConds(cond)));
        thisVisStimOrder = RPiVisAllMap.Var2(matches(RPiVisAllMap.Var1,visMapConds(cond)));
        % This is a N by 2 matrix, where first row holds the identity
        % of the grating and the second row holds the timestamps
        sortedThisVisStimOrder = sortrows([thisVisStimOrder thisVisStimTimes]); 

        assert(length(sortedThisVisStimOrder) == numOfStims) % just make sure extraction was correct
        
        temp_fullDuration = [];
        temp_partialDuration = [];
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
            %numOfSpikesOrientation(cond,thisVisStimIdentity,:) = numOfSpikesOrientation(cond,thisVisStimIdentity,:) + N;
            
            % for partial
            indivTrial = alignSpikes(alignSpikes <= 1.5 & alignSpikes >= 0.5);
            N = histcounts(indivTrial, edgesPartial);
            temp_partialDuration = [temp_partialDuration; N];
            %numOfSpikesOrientation_partial(cond,thisVisStimIdentity,:) = numOfSpikesOrientation_partial(cond,thisVisStimIdentity,:) + N;
        end

        numOfSpikesOrientation{cond,1} = temp_fullDuration;
        numOfSpikesOrientation_partial{cond,1} = temp_partialDuration;
        numOfSpikesOrientation{cond,2} = sortedThisVisStimOrder(:,1);
        numOfSpikesOrientation_partial{cond,2} = sortedThisVisStimOrder(:,1);

        % Plot other aspects of the scatter (line at 0, etc)
        plot(zeros(1, numOfStims+1), 1:numOfStims+1, 'k--', 'LineWidth', 1.0)
        hold on;
        fill([0, 1.5, 1.5, 0], [0, 0, numOfStims, numOfStims], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
        xlim([-0.3, 2.0]);
        ylim([0, numOfStims+1]);
        xlabel('Time (s)');
        ylabel('Ordered by orientation');
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10;
    end

    
    %% Last, plotting the orientation plots 
    % numOfSpikesOrientation is a cell array, where the first column
    % holds the PSTH for each stimulus presentation
    % second column holds the stimulus order matrix


    sumSpikesOrientation = cellfun(@(x) sum(x,2), numOfSpikesOrientation, 'UniformOutput', false);
    sumSpikesOrientation_partial = cellfun(@(x) sum(x,2), numOfSpikesOrientation_partial, 'UniformOutput', false);
    statsOrientation = {}; % this is a cell array where each row is a condition and column is a visual map identity. 
    statsOrientation_partial = {};
    % First element of the matrix is the mean and the second element is
    % std, 3rd element is SEM
   
    for cond = 1:length(visMapConds)
        for j = 1:length(visMapIdentity)
            temp_ind = find(sumSpikesOrientation{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation{cond,1}(temp_ind)/(max(edgesOrientation)-min(edgesOrientation));
            statsOrientation{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];
        
            temp_ind = find(sumSpikesOrientation_partial{cond,2} == visMapIdentity(j));
            thisIdentitySpikes = sumSpikesOrientation_partial{cond,1}(temp_ind)/(max(edgesOrientation)-min(edgesOrientation));
            statsOrientation_partial{cond,j} = [mean(thisIdentitySpikes) std(thisIdentitySpikes) std(thisIdentitySpikes)/sqrt(length(thisIdentitySpikes))];     
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
        
        fullSEM = cellfun(@(x) x(3), statsOrientation);
        fullMean = cellfun(@(x) x(1), statsOrientation);
        partialSEM = cellfun(@(x) x(3), statsOrientation_partial);
        partialMean = cellfun(@(x) x(1), statsOrientation_partial);

        maxYlim = max(fullMean,[],'all');

        yyaxis left
        errorbar(1:30:360, fullMean(cond,:), fullSEM(cond,:),'b-','LineWidth',1.4,'Capsize',0);
        ylabel('Full 1.5 seconds')
        ylim([0 maxYlim+3])
        hold on;
        yyaxis right
        errorbar(1:30:360, partialMean(cond,:), partialSEM(cond,:),'r-','LineWidth',1.4,'Capsize',0);
        ylabel('Excluding first 0.5 seconds')
        title(append('Response to ',visMapConds{cond}));
        ax = gca;
        ax.FontSize = 10; 
        xlim([0 360])
        ylim([0 maxYlim+3])
    end

    sgtitle(['Cluster ',num2str(goodClusters(i)+1)])
    fPath = fullfile(myKsDir,'Orientation');
    try
        saveFigurePdf(unitSpikePSTHScatter, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
    catch ME
        disp(['failed to save figure'])
    end
    close(unitSpikePSTHScatter);

end

end