%% Set up paths to import the data
format longg 

myKsDir = cd; % cd to the kilosort output folder
Ksparams.excludeNoise = false; % we make sure all the spike times are imported. 
sp = loadKSdir(myKsDir,Ksparams); % we load the structure in.
spiketimes = readmatrix('spike_times.txt');  % this is converted using the accurate spike rate time
sp.st = spiketimes; 
sp.ss = readNPY('spike_times.npy'); % in samples


goodClu_idx = find(sp.cgs == 2); % find index of good clusters
goodClusters = sp.cids(goodClu_idx); % this is -1 from the MATLAB indexing

% Load visual stimulation timestamps
bindataDir = fullfile(cd, '..','*_1500_TPrime.txt');
visStimFiles = dir(bindataDir);
visStimTimes = readmatrix(fullfile(visStimFiles.folder,visStimFiles.name));

% Load the order of visual stimulation
eventDir= fullfile(cd, '..','*_visStimOrder.txt');
eventDirFiles = dir(eventDir);
eventTable = readtable(fullfile(eventDirFiles.folder,eventDirFiles.name),'Delimiter',{':','-'});
 
RPiVisAllMap = eventTable(contains(eventTable.Var1,'map'),:);
RPiVisAllMap.StimTimes = visStimTimes(1:height(RPiVisAllMap));

RPiVisSineMap = RPiVisAllMap(matches(RPiVisAllMap.Var1,'Sine map'),:);
RPiVisGratingMap = RPiVisAllMap(matches(RPiVisAllMap.Var1,'Grating map'),:);
RPiVisSquareMap = RPiVisAllMap(matches(RPiVisAllMap.Var1,'Square map'),:);

%% Set up an array that describes the visual stimulation protocol
% Protocol is usually 1.5 seconds of visual stimulus and 2 seconds of ISI
assert(mod(length(visStimTimes),120) == 0) % Check that visual stimulation time is a multiple of 120

edges = -0.3:0.005:2.0;
numOfbins = length(edges)-1;

plotVisStimScatter(sp, goodClusters, visStimTimes, myKsDir, RPiVisAllMap)
close all;
plotVisStimOptoScatter(sp, goodClusters, visStimTimes, myKsDir, RPiVisAllMap)
close all;
calculateFRComparison(sp, goodClusters, myKsDir, RPiVisAllMap)


% 
% for i = 1:length(goodClusters)
%     unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds
%     numOfSpikesPSTH = zeros(1,numOfbins);
% 
%     unitSpikePSTHScatter = figure('Renderer', 'painters', 'Position', [10 10 1200 1200]);
%     subplot(2,2,1)
%     y = 0;
%     for stim = 1:120
%         scatter(unitSpikes - visStimTimes(stim),ones(1,length(unitSpikes)) + y, '.b')        % Align the spike train to the pulse
%         y = y + 1;                                                                          % Stack them on top of each other
%         hold on
% 
%         alignSpikes = unitSpikes - visStimTimes(stim);  
%         indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3); % Play with the lower minima to obtain baseline
%         N = histcounts(indivTrial,edges);
%         numOfSpikesPSTH = numOfSpikesPSTH + N;
%     end
%     plot(zeros(1, 120), 1:120, 'k--', 'LineWidth', 1.0)     % Plot a dashed line at the onset of the spike
%     hold on;
%     fill([0 1.5 1.5 0],[0 0 length(visStimTimes) length(visStimTimes)],'r','FaceAlpha',0.1,'EdgeColor','none')
%     xlim([-0.3,2.0])
%     ylim([0 121])
%     xlabel('Time (s)')
%     ylabel('Trial')
%     title('Response to sine grating')
%     ax = gca;
%     ax.FontSize = 10;
% 
%     subplot(2,2,2)
%     spikeRatePSTH = numOfSpikesPSTH/(240*0.005);
%     maxYlim = max(spikeRatePSTH);
%     plot(spikeRatePSTH)
%     xticks(1:50:451)
%     xticklabels(-0.3:0.25:2.0)
%     hold on;
%     fill([61 361 361 61],[0 0 maxYlim+3 maxYlim+3],'r','FaceAlpha',0.1,'EdgeColor','none')
%     ylim([0 maxYlim+3])
%     xlabel('Time (s)')
%     ylabel('Hz')
%     title('Response to sine grating')
% 
%     numOfSpikesPSTH = zeros(1,numOfbins);
%     % Square grating now
%     subplot(2,2,3)
%     y = 0;
%     for stim = 121:240
%         scatter(unitSpikes - visStimTimes(stim),ones(1,length(unitSpikes)) + y, '.b')        % Align the spike train to the pulse
%         y = y + 1;                                                                          % Stack them on top of each other
%         hold on
% 
%         alignSpikes = unitSpikes - visStimTimes(stim);  
%         indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3); % Play with the lower minima to obtain baseline
%         N = histcounts(indivTrial,edges);
%         numOfSpikesPSTH = numOfSpikesPSTH + N;
%     end
%     plot(zeros(1, 240), 1:240, 'k--', 'LineWidth', 1.0)     % Plot a dashed line at the onset of the spike
%     hold on;
%     fill([0 1.5 1.5 0],[0 0 length(visStimTimes) length(visStimTimes)],'r','FaceAlpha',0.1,'EdgeColor','none')
%     xlim([-0.3,2.0])
%     ylim([0 121])
%     xlabel('Time (s)')
%     ylabel('Trial')
%     title('Response to square grating')
%     ax = gca;
%     ax.FontSize = 10;
% 
%     subplot(2,2,4)
%     spikeRatePSTH = numOfSpikesPSTH/(240*0.005);
%     maxYlim = max(spikeRatePSTH);
%     plot(spikeRatePSTH)
%     xticks(1:50:451)
%     xticklabels(-0.3:0.25:2.0)
%     hold on;
%     fill([61 361 361 61],[0 0 maxYlim+3 maxYlim+3],'r','FaceAlpha',0.1,'EdgeColor','none')
%     ylim([0 maxYlim+3])
%     xlabel('Time (s)')
%     ylabel('Hz')
%     title('Response to square grating')
% 
%     sgtitle(['Cluster ',num2str(goodClusters(i)+1)])
% 
%     fPath = fullfile(myKsDir,'scatter');
%     saveFigurePdf(unitSpikePSTHScatter, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
%     close all
% end
% 
% 
% %% Now I want to calculate orientation selectivity index
% 
% % Lien and Scanziani, 2018 Nature (https://www.nature.com/articles/s41586-018-0148-5#Sec7)
% % - they used 2.3s long grating and calculated the OSI using 2.0 seconds
% % (removed the first 0.3s from the stimulus presentation)
% % - used average spike rate during the response window
% %
% % Velez-Fort et al., 2014 Neuron
% % (https://www.sciencedirect.com/science/article/pii/S089662731400676X?via%3Dihub#abs0015)
% % - Analyzed the second half of the 2 second drift (this removes
% % non-specific direction onset responses)
% 
% %%%%%%%%%%% to do: make the indices, number of stim presentations, all NOT
% %%%%%%%%%%% hard coded
% 
% SineOrder = [RPiVisSineMap.Var2 RPiVisSineMap.StimTimes];
% sortedSineOrder = sortrows(SineOrder);
% SquareOrder = [RPiVisSquareMap.Var2 RPiVisSquareMap.StimTimes];
% sortedSquareOrder = sortrows(SquareOrder);
% GratingOrder = [RPiVisGratingMap.Var2 RPiVisGratingMap.StimTimes];
% sortedGratingOrder = sortrows(GratingOrder);
% 

% for i = 1:length(goodClusters)
%     unitSpikes = sp.st(sp.clu == goodClusters(i)); % already in seconds
%     numOfSpikesPSTH = zeros(1,numOfbins);
% 
%     unitSpikePSTHScatter = figure('Renderer', 'painters', 'Position', [10 10 600 1200]);
%     subplot(2,1,1)
%     y = 0;
%     for stim = 1:length(sortedSineOrder)
%         alignSpikes = unitSpikes - sortedSineOrder(stim,2);  
%         scatter(alignSpikes,ones(1,length(unitSpikes)) + y, '.b')        % Align the spike train to the pulse
%         y = y + 1;                                                                          % Stack them on top of each other
%         hold on
%         indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3); % Play with the lower minima to obtain baseline
%         N = histcounts(indivTrial,edges);
%         numOfSpikesPSTH = numOfSpikesPSTH + N;
%     end
%     plot(zeros(1, 120), 1:120, 'k--', 'LineWidth', 1.0)     % Plot a dashed line at the onset of the spike
%     hold on;
%     fill([0 1.5 1.5 0],[0 0 length(sortedSineOrder) length(sortedSineOrder)],'r','FaceAlpha',0.1,'EdgeColor','none')
%     xlim([-0.3,2.0])
%     ylim([0 121])
%     xlabel('Time (s)')
%     ylabel('Trial')
%     title('Response to sine grating')
%     ax = gca;
%     ax.FontSize = 10;
% 
%     subplot(2,1,2)
%     y = 0;
%     for stim = 1:length(sortedGratingOrder)
%         alignSpikes = unitSpikes - sortedGratingOrder(stim,2);  
%         scatter(alignSpikes,ones(1,length(unitSpikes)) + y, '.b')        % Align the spike train to the pulse
%         y = y + 1;                                                                          % Stack them on top of each other
%         hold on
%         indivTrial = alignSpikes(alignSpikes < 2.0 & alignSpikes >= -0.3); % Play with the lower minima to obtain baseline
%         N = histcounts(indivTrial,edges);
%         numOfSpikesPSTH = numOfSpikesPSTH + N;
%     end
%     plot(zeros(1, 120), 1:120, 'k--', 'LineWidth', 1.0)     % Plot a dashed line at the onset of the spike
%     hold on;
%     fill([0 1.5 1.5 0],[0 0 length(sortedGratingOrder) length(sortedGratingOrder)],'r','FaceAlpha',0.1,'EdgeColor','none')
%     xlim([-0.3,2.0])
%     ylim([0 121])
%     xlabel('Time (s)')
%     ylabel('Trial')
%     title('Response to square grating')
%     ax = gca;
%     ax.FontSize = 10;
% 
%     sgtitle(['Cluster ',num2str(goodClusters(i)+1)])
%     fPath = fullfile(myKsDir,'Orientation');
%     saveFigurePdf(unitSpikePSTHScatter, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
%     close all
% end



