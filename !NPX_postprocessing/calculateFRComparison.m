function calculateFRComparison(sp, goodClusters, myKsDir, RPiVisAllMap)
% calculateFRComparison(sp, goodClusters, myKSDir, RPiVisAllMap)
% sp: structure that holds ks outputs
% goodClusters: good cluster assignment
% visStimTimes: visual stimulation timestamps (.txt file, usually in a
% table)
% myKsDir: pathway to kilosort
% RPiVisAllMap: table that holds visual directions, as well as the
% timestamps
%
% This code extracts the baseline firing rate (1.5s prior to visual
% stimulation onset) vs. the active firing rate (1.5s of visual
% presentation) 
% It also plots the results in a scatter (with y=x unity line).
% 
%
% Written on 8/19/2024 
% Chris (Seong Yeol) An
% san11@jh.edu

%% We first adjust the RPIVisAllMap table to identify the visual stimulation epochs with no interleaved opto stimulations.
% These epochs will be marked with ' - 1 '
visMapConds = unique(RPiVisAllMap.Var1); % These are list of visual mappings that the animal viewed in this session
visMapIdentity = unique(RPiVisAllMap.Var2);
pureVisMapStimTimes = [];
for i = 1:length(RPiVisAllMap.Var3)
    temp = RPiVisAllMap.Var3(i);
    if isnan(temp)
        RPiVisAllMap.Var3(i) = -1; % replace the empty elements with ' - 1 '
        pureVisMapStimTimes = [pureVisMapStimTimes RPiVisAllMap.StimTimes(i)]; % we are extracting timestamps
    end
end

%% We will not iterate through the good clusters and extract corresponding active and baseline FRs
% First, initialize matrices to hold active and baseline FR.
activeFR = zeros(length(goodClusters),3);
baselineFR = zeros(length(goodClusters),3);
% We make the first column of this matrix the cluster number (for future
% identification). The second column will hold the mean, and 3rd column
% SEM
activeFR(:,1) = goodClusters + 1;
baselineFR(:,1) = goodClusters + 1;
% Now we will iterate through goodClusters
for i = 1:length(goodClusters)
    % First extract the corresponding unit spikes
    unitSpikes = sp.st(sp.clu == goodClusters(i)); % in seconds
    % Initialize temporary matrices to hold the Hz
    baseline_clu = zeros(length(pureVisMapStimTimes),1);
    active_clu = zeros(length(pureVisMapStimTimes),1);
    % Then we will align the spikes while iterating through the spiketimes
    for stim = 1:length(pureVisMapStimTimes)
        alignedSpikes = unitSpikes - pureVisMapStimTimes(stim); % align the spikes to stimulation onset
        relevantSpikes = alignedSpikes(alignedSpikes >= -1.5 & alignedSpikes <= 1.5); % extract the spikes in baseline vs. active
   
        % Now I want to calculate baseline and active FR and add to the
        % temporary matrix
        baseline_clu(stim) = length(find(relevantSpikes <= 0))/1.5;
        active_clu(stim) = length(find(relevantSpikes > 0))/1.5;  
    end
   
    % Now I want to calculate mean and SEM
    baselineFR(i,2) = mean(baseline_clu);
    baselineFR(i,3) = std(baseline_clu)/sqrt(length(baseline_clu));
    activeFR(i,2) = mean(active_clu);
    activeFR(i,3) = std(active_clu)/sqrt(length(active_clu));
end
% I want to find which points are above the unity line.
boolEngagedByVisStim = zeros(length(goodClusters),1);
for clu = 1:length(goodClusters)
    if activeFR(clu,2) > baselineFR(clu,2)
        boolEngagedByVisStim(clu) = 1;
    else
        boolEngagedByVisStim(clu) = 0;
    end
end
percEngagedByVisStim = sum(boolEngagedByVisStim)/length(boolEngagedByVisStim) * 100;
% Calculate VisStimModulationIndex
% We calculate this by (active-baseline)/active+baseline
VMI = (activeFR(:,2) - baselineFR(:,2))./(activeFR(:,2) + baselineFR(:,2));

% Convert matrices into tables
activeFR_Table = array2table(activeFR,'VariableNames',{'Cluster','Mean','SEM'});
baselineFR_Table = array2table(baselineFR,'VariableNames',{'Cluster','Mean','SEM'});


%% Now draw the relevant plots
activeBaselineFRPlot = figure('Renderer', 'painters', 'Position', [10 10 1200 1200]);
% Plot scatter
subplot(2,6,[2 3])
scatter(baselineFR(:,2), activeFR(:,2),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
hold on;
plot(0:50,0:50,'-r','LineWidth',1.5)
xlim([0 50]);
ylim([0 50]);
xticks(0:10:50);
yticks(0:10:50);
xlabel('Baseline (Hz)')
ylabel('Active (Hz)')
axis square
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1.3;
title(append(num2str(percEngagedByVisStim),'%'))
% Plot y-histogram
subplot(2,8,1)
binEdges = 0:1:50;
histogram(activeFR(:,2),binEdges,'orientation','horizontal','Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[0 0 0])
set(gca,'xdir','reverse')
box off;
axis off;
% Plot x-histogram
subplot(4,6,[14 15])
histogram(baselineFR(:,2),binEdges,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[0 0 0])
set(gca,'ydir','reverse')
box off;
axis off;

subplot(2,6,[11 12])
binEdges = -1.5:0.1:1.5;
histogram(VMI,binEdges,'Normalization','percentage','FaceColor',[0 0 0],'EdgeColor',[0 0 0])
hold on;
xline(nanmean(VMI),'r-','LineWidth',2.0)
xlabel('Visual Modulation Index')
ylabel('Percentage')
title(append('mean = ',num2str(nanmean(VMI))))
axis square
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1.3;
box off

%% Save the file
parentFolderPath = fileparts(myKsDir);
[~,parentFolderName] = fileparts(parentFolderPath);
sgtitle(parentFolderName)
fPath = fullfile(myKsDir,'Analysis_plots',append(parentFolderName,'_BaselineActiveFR'));
saveFigurePdf(activeBaselineFRPlot,fPath)
pause(10);
close;

% Save the variables
fname = fullfile(myKsDir,'analysis_variables','VMI_firingRates');
save(fname,'activeFR_Table','baselineFR_Table','VMI');



end