function calculateOSI(goodClusters, myKsDir)
% myKsDir: pathway to kilosort
% goodClusters: good cluster assignment
% This code calculates OSI for every cluster
%
% Written on 8/19/2024 
% Chris (Seong Yeol) An
% san11@jh.edu

%% We should iterate through the clusters and calculate OSI
parentPath = fullfile(myKsDir,'OptoOrientation','struct');
OSI = [];
OSI_margrie = [];
DSI_margrie = [];
for i = 1:length(goodClusters)
    cluName = append('Cluster_',num2str(goodClusters(i)+1),'.mat');
    loadPath = fullfile(parentPath,cluName);
    load(loadPath); % this loads the structure orientationStats

    % structure orientationStats
    % cluster: cluster number
    % fullMean: mean values for each visual direction
    % fullSEM: same as above but SEM
    % numSpikes: raw values that were used to calculate the mean
    % binEdges: edges for the bins in numSpikes
    % visDirections: visual mapping directions (usually 12)

    %% To calculate OSI, we use the following formula
    % OSI = sqrt((sum(r_k*sin(2h_k))^2+(sum(r_k*cos(2h_k)))^2)/sum(r_k), where
    % r_k is the response to the kth direction by direction h_k
    % Olsen et al., 2012 Nature; Kerlin et al., 2010 Neuron; 
    % Ringach et al.,2022 J. Neurosci.
    
    % First, we extract h_k. This is held in visDirections
    h_k = orientationStats.visDirections;
    h_k_rad = deg2rad(h_k);
    % Then, we will now calculate r_k, which is the average response to h_k
    visConds = fieldnames(orientationStats.numSpikes); % Visual conditions (usually pure, opto, and non-opto)
    r_k = zeros(length(visConds),length(h_k)); % initialize r_k for average response
    % numSpikes has elements that are equal to the binEdges
    % to calculate r_k, we need to find the indices where the visual
    % stimulation onset happens. This is from timestamp = 0 to timestamp = 0 in
    % the edges
    startTimeIdx = find(orientationStats.binEdges <= 0.50001 & orientationStats.binEdges >= 0.49999);
    endTimeIdx = find(orientationStats.binEdges <= 1.00001 & orientationStats.binEdges >= 0.99999)-1;
    % Iterate through the conditions
    for conds = 1:length(visConds)
        tmp = orientationStats.numSpikes.(visConds{conds}){1,1};
        tempVisIdentity = orientationStats.numSpikes.(visConds{conds}){1,2}; % visMapIdentity
        tempSpikeCounts = tmp(:,startTimeIdx:endTimeIdx); % this is a matrix size of tempVisIdentity x binEdges
        % Now I want to calculate r_k for each visual direction
        % We iterate through each direction...
        for visDis = 1:length(h_k)
            tempDirectionSpikeCounts = tempSpikeCounts(tempVisIdentity == h_k(visDis),:);
            tempDirectionSpikeRate = sum(tempDirectionSpikeCounts,2)/0.5; % Spike rate for this given direction
            r_k(conds,visDis) = mean(tempDirectionSpikeRate);
        end
    end

    % Now use the formula to calculate OSI for each condition
    OSI_clu = zeros(length(visConds),1);
    for conds = 1:length(visConds)
        sinV = sum(r_k(conds,:)'.*sin(2*h_k_rad));
        cosV = sum(r_k(conds,:)'.*cos(2*h_k_rad));
        numerator = sqrt(sinV^2 + cosV^2);
        denominator = sum(r_k(conds,:));
        OSI_clu(conds,1) = numerator/denominator;
    end
    % Add to OSI matrix
    OSI = [OSI; OSI_clu'];

    %% Try another OSI approach
    % Find the vector average of the responses for the 12 directions, and
    % use the direction closest to the value of the vector average as the
    % preferred direction.
    % Velez-Fort et al., 2014
    % First convert the angles into rad
    h_k_rad = h_k' * pi / 180;
    mean_mag = zeros(length(visConds),1);
    mean_direction_deg = zeros(length(visConds),1);
    mean_direction_rad = zeros(length(visConds),1);
    OSI_margrie_clu = zeros(length(visConds),1);
    DSI_margrie_clu = zeros(length(visConds),1);
    for conds = 1:length(visConds)
        % Calculate the x and the y elements of the vector.
        x_vector = r_k(conds,:) .* cos(h_k_rad);
        y_vector = r_k(conds,:) .* sin(h_k_rad);
        % Sum each vector component
        x_vectorSum = sum(x_vector);
        y_vectorSum = sum(y_vector);
        % the mean magnitude and the direction is the resulting magnitude and
        % direction
        mean_mag(conds,1) = hypot(x_vectorSum,y_vectorSum);
        mean_direction_rad(conds,1) = atan2(y_vectorSum,x_vectorSum);
        mean_direction_deg(conds,1) = (atan2(y_vectorSum,x_vectorSum)) * 180/pi;
        if (atan2(y_vectorSum,x_vectorSum)) * 180/pi < 0
            mean_direction_deg(conds,1) = (atan2(y_vectorSum,x_vectorSum)) * 180/pi + 360;
        end
        % Now find the closest angle to the mean_direction_deg
        [~,tmpIdx] = min(abs(h_k - (atan2(y_vectorSum,x_vectorSum)) * 180/pi));
        preferredOrientation = h_k(tmpIdx);
        orthoOrientation = preferredOrientation - 90;
        nullOrientation = preferredOrientation + 180;
        if orthoOrientation < 0
            orthoOrientation = orthoOrientation + 360;
        end
        if nullOrientation >= 360
            nullOrientation = nullOrientation - 360;
        end
        % OSI = (pref-ortho)/(pref+ortho); DSI = (pref-null)/(pref+null)
        orthoOrientationIdx = find(h_k == orthoOrientation);
        preferredOrientationIdx = find(h_k == preferredOrientation);
        nullOrientationIdx = find(h_k == nullOrientation);
        OSI_margrie_clu(conds,1) = (r_k(conds,preferredOrientationIdx) - r_k(conds,orthoOrientationIdx))./...
            (r_k(conds,preferredOrientationIdx) + r_k(conds,orthoOrientationIdx));
        DSI_margrie_clu(conds,1) = (r_k(conds,preferredOrientationIdx) - r_k(conds,nullOrientationIdx))./...
            (r_k(conds,nullOrientationIdx) + r_k(conds,nullOrientationIdx));
    end
    % Add to the matrix
    OSI_margrie = [OSI_margrie; OSI_margrie_clu'];
    DSI_margrie = [DSI_margrie; DSI_margrie_clu'];

    %% Generate polar plot figure (like those from Velez-Fort et al., 2014)
    figColors = [1 0 0;0 0 1;0 0 0;0 1 0];
    polarOrientationFig = figure('Renderer', 'painters', 'Position', [10 10 600 600]);
    legTxt = visConds;
    OSItxt = visConds;
    for conds = 1:length(visConds)
        polarplot([h_k_rad h_k_rad(1)],[r_k(conds,:) r_k(conds,1)],'Color',figColors(conds,:),'LineWidth',2)
        hold on;
        polarplot([0 mean_direction_rad(conds)], [0 mean_mag(conds)],'Color',figColors(conds,:), 'LineWidth', 3.0, 'HandleVisibility','off');
        hold on;
        legTxt{conds} = append(legTxt{conds},"=",num2str(OSI_clu(conds)));
        OSItxt{conds} = append(OSItxt{conds},"=",num2str(OSI_margrie_clu(conds)));
    end
    annotation('textbox',[0.7, 0.2, 0, 0],'string',OSItxt)
    legend(legTxt)
    sgtitle(['Cluster ',num2str(goodClusters(i)+1)])
    ax = gca;
    ax.FontSize = 10;
    ax.LineWidth = 1.3;

    %% Save polar plot figure
    fPath = fullfile(myKsDir,'PolarOrientation');
    saveFigurePdf(polarOrientationFig, fullfile(fPath, append('Cluster_',num2str(goodClusters(i)+1))))
end
close all;
%% Save OSI results 
OSI_Table = array2table(OSI,'VariableNames',visConds);
rowNames = arrayfun(@num2str,goodClusters+1,'UniformOutput',false);
OSI_Table.Cluster = rowNames';
fname = fullfile(myKsDir,'analysis_variables','OSI_Table');
save(fname,'OSI_Table')

OSI_margrie_Table = array2table(OSI_margrie,'VariableNames',visConds);
rowNames = arrayfun(@num2str,goodClusters+1,'UniformOutput',false);
OSI_margrie_Table.Cluster = rowNames';
fname = fullfile(myKsDir,'analysis_variables','OSI_margrie_Table');
save(fname,'OSI_margrie_Table')

%% plot scatters
% % first find the proper index for opto and no-opto
% optoIdx = find(matches(visConds,"opto") == 1);
% noOptoIdx = find(matches(visConds,"noOpto") == 1);
% 
% OSIScatter = figure('Renderer', 'painters', 'Position', [10 10 1200 600]);
% subplot(1,2,1)
% scatter(OSI(:,noOptoIdx),OSI(:,optoIdx),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
% hold on;
% plot(0:0.1:1,0:0.1:1,'r-')
% axis square
% box off
% xlabel('No-Opto')
% ylabel('Opto')
% title('OSI from Olsen 2012')
% 
% subplot(1,2,2)
% scatter(OSI_margrie(:,noOptoIdx),OSI_margrie(:,optoIdx),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
% hold on;
% plot(0:0.1:1,0:0.1:1,'r-')
% axis square
% box off
% xlabel('No-Opto')
% ylabel('Opto')
% title('OSI from Velez-Fort 2014')


end