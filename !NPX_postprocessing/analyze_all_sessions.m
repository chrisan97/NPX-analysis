%% Load OSI
% Initialize an empty structure to store the loaded data
OSI_struct = struct();

% List immediate subfolders of rootDir
rootDir = cd;
subfolders = dir(rootDir);
subfolderNames = {subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'})).name};
targetFileName = 'OSI_Table'

% Start the recursive search from the root directory
OSI_struct = searchAndLoad(rootDir, targetFileName, OSI_struct, rootDir);

%% Load waveform width
wfWidth_struct = struct();

% List immediate subfolders of rootDir
rootDir = cd;
subfolders = dir(rootDir);
subfolderNames = {subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'})).name};
targetFileName = 'waveformWidth'

% Start the recursive search from the root directory
wfWidth_struct = searchAndLoad(rootDir, targetFileName, wfWidth_struct, rootDir);

%% Load waveform depth
wfDepth_struct = struct();

% List immediate subfolders of rootDir
rootDir = cd;
subfolders = dir(rootDir);
subfolderNames = {subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'})).name};
targetFileName = 'waveformDepth'

% Start the recursive search from the root directory
wfDepth_struct = searchAndLoad(rootDir, targetFileName, wfDepth_struct, rootDir);

%% Load optoModulation
OMI_struct = struct();

% List immediate subfolders of rootDir
rootDir = cd;
subfolders = dir(rootDir);
subfolderNames = {subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'})).name};
targetFileName = 'optoModulationTable'

% Start the recursive search from the root directory
OMI_struct = searchAndLoad(rootDir, targetFileName, OMI_struct, rootDir);

%% Load visualModulation
VMI_struct = struct();

% List immediate subfolders of rootDir
rootDir = cd;
subfolders = dir(rootDir);
subfolderNames = {subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'})).name};
targetFileName = 'VMI_firingRates'

% Start the recursive search from the root directory
VMI_struct = searchAndLoad(rootDir, targetFileName, VMI_struct, rootDir);


%% save the combined structure
save('090224-combined_struct','OMI_struct','OSI_struct','wfDepth_struct','wfWidth_struct','VMI_struct')


%% First plot the OSI of noOpto to pure
fields = fieldnames(OSI_struct);

combinedNoOpto = [];
combinedOpto = [];
combinedPure = [];
for i = 1:length(fields)
    fieldName = fields{i};
    thisTable = OSI_struct.(fieldName).OSI_Table;

    combinedNoOpto = [combinedNoOpto; thisTable.noOpto];
    combinedOpto = [combinedOpto; thisTable.opto];
    combinedPure = [combinedPure; thisTable.pure];

    thisOSIPlots = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);
    subplot(3,2,1)
    scatter(thisTable.pure,thisTable.noOpto,10,'k','filled')
    hold on
    plot(0:0.1:1,0:0.1:1,'r-')
    xlabel('Pure vis stim')
    ylabel('Interleaved no-Opto stim')
    xlim([0 1])
    ylim([0 1])
    axis square
    ax = gca;
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    title('Pure vs noOpto')
    box off;

    subplot(3,2,2)
    mdl = fitlm(thisTable.pure,thisTable.noOpto);
    hold on;
    plot(mdl)
    xlabel('Pure vis stim')
    ylabel('Interleaved no-Opto stim')
    xlim([0 1])
    ylim([0 1])
    axis square
    legend off;
    ax = gca;
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    box off;

    subplot(3,2,3)
    scatter(thisTable.pure,thisTable.opto,10,'k','filled')
    hold on
    plot(0:0.1:1,0:0.1:1,'r-')
    xlabel('Pure vis stim')
    ylabel('Interleaved Opto stim')
    xlim([0 1])
    ylim([0 1])
    axis square
    ax = gca;
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    title('Pure vs Opto')
    box off;

    subplot(3,2,4)
    mdl = fitlm(thisTable.pure,thisTable.opto);
    hold on;
    plot(mdl)
    xlabel('Pure vis stim')
    ylabel('Interleaved Opto stim')
    xlim([0 1])
    ylim([0 1])
    axis square
    legend off;
    ax = gca;
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    box off;

    subplot(3,2,5)
    scatter(thisTable.noOpto,thisTable.opto,10,'k','filled')
    hold on
    plot(0:0.1:1,0:0.1:1,'r-')
    xlabel('Interleaved no-Opto stim')
    ylabel('Interleaved Opto stim')
    xlim([0 1])
    ylim([0 1])
    axis square
    legend off;
    ax = gca;
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    title('noOpto vs Opto')
    box off;

    subplot(3,2,6)
    mdl = fitlm(thisTable.noOpto,thisTable.opto);
    hold on;
    plot(mdl)
    xlabel('Interleaved no-Opto stim')
    ylabel('Interleaved Opto stim')
    xlim([0 1])
    ylim([0 1])
    axis square
    legend off;
    ax = gca;
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    box off;

    sgtitle(fieldName)
    fPath = fullfile(cd,'combined_analysis');
    try
        saveFigurePdf(thisOSIPlots, fullfile(fPath, append(fieldName,'_OSI')))    
    catch ME
        disp(['failed to save figure'])
    end
end

% Combined
combinedOSIScatter = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);
subplot(3,2,1)
scatter(combinedPure,combinedNoOpto,10,'k','filled')
hold on
plot(0:0.1:1,0:0.1:1,'r-')
xlabel('Pure vis stim')
ylabel('Interleaved no-Opto stim')
xlim([0 1])
ylim([0 1])
axis square
ax = gca;
xticks(0:0.2:1)
yticks(0:0.2:1)
ax.FontSize = 14;
title('Pure vs noOpto')
box off;

subplot(3,2,2)
mdl = fitlm(combinedPure,combinedNoOpto);
hold on;
plot(mdl)
xlabel('Pure vis stim')
ylabel('Interleaved no-Opto stim')
xlim([0 1])
ylim([0 1])
axis square
legend off;
ax = gca;
xticks(0:0.2:1)
yticks(0:0.2:1)
ax.FontSize = 14;
box off;

subplot(3,2,3)
scatter(combinedPure,combinedOpto,10,'k','filled')
hold on
plot(0:0.1:1,0:0.1:1,'r-')
xlabel('Pure vis stim')
ylabel('Interleaved Opto stim')
xlim([0 1])
ylim([0 1])
axis square
ax = gca;
xticks(0:0.2:1)
yticks(0:0.2:1)
ax.FontSize = 14;
title('Pure vs Opto')
box off;

subplot(3,2,4)
mdl = fitlm(combinedPure,combinedOpto);
hold on;
plot(mdl)
xlabel('Pure vis stim')
ylabel('Interleaved Opto stim')
xlim([0 1])
ylim([0 1])
axis square
legend off;
ax = gca;
xticks(0:0.2:1)
yticks(0:0.2:1)
ax.FontSize = 14;
box off;

subplot(3,2,5)
scatter(combinedNoOpto,combinedOpto,10,'k','filled')
hold on
plot(0:0.1:1,0:0.1:1,'r-')
xlabel('Interleaved no-Opto stim')
ylabel('Interleaved Opto stim')
xlim([0 1])
ylim([0 1])
axis square
legend off;
ax = gca;
xticks(0:0.2:1)
yticks(0:0.2:1)
ax.FontSize = 14;
title('noOpto vs Opto')
box off;

subplot(3,2,6)
mdl = fitlm(combinedNoOpto,combinedOpto);
hold on;
plot(mdl)
xlabel('Interleaved no-Opto stim')
ylabel('Interleaved Opto stim')
xlim([0 1])
ylim([0 1])
axis square
legend off;
ax = gca;
xticks(0:0.2:1)
yticks(0:0.2:1)
ax.FontSize = 14;
box off;

sgtitle('Combined OSI')
saveFigurePdf(combinedOSIScatter, fullfile(fPath, append(fieldName,'_combinedOSI')))    
close all;

%% Just waveform
fields = fieldnames(wfWidth_struct);
combinedwfWidth = [];
for i = 1:length(fields)
    fieldName = fields{i};
    thiswfWidth = wfWidth_struct.(fieldName).waveformWidth(2,:);
    combinedwfWidth = [combinedwfWidth thiswfWidth];

    thisWfWidthPlots = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
    histogram(thiswfWidth,0:0.1:2.0,'FaceColor','black');
    xlabel('Waveform width (ms)')
    ylabel('Counts')
    ax = gca;
    ax.FontSize = 14;
    box off;

    sgtitle(fieldName)
    fPath = fullfile(cd,'combined_analysis');
    try
        saveFigurePdf(thisWfWidthPlots, fullfile(fPath, append(fieldName,'_wfWidth')))    
    catch ME
        disp(['failed to save figure'])
    end
end

combinedWfWidthPlots = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
histogram(combinedwfWidth,0:0.035:2.0,'FaceColor','black');
xlabel('Waveform width (ms)')
ylabel('Counts')
box off;
set(gca,'LineWidth',2,'FontSize',14)
sgtitle('Combined waveform width')
saveFigurePdf(combinedWfWidthPlots, fullfile(fPath, append(fieldName,'_combinedWfWidth')))    
close all;

%% OSI and waveform
fields = fieldnames(wfWidth_struct);

combinedwfWidth = [];
combinedNoOpto = [];
combinedOpto = [];
combinedPure = [];

for i = 1:length(fields)
    fieldName = fields{i};
    thisTable = OSI_struct.(fieldName).OSI_Table;
    thiswfWidth = wfWidth_struct.(fieldName).waveformWidth(2,:);
    combinedwfWidth = [combinedwfWidth thiswfWidth];
    combinedNoOpto = [combinedNoOpto; thisTable.noOpto];
    combinedOpto = [combinedOpto; thisTable.opto];
    combinedPure = [combinedPure; thisTable.pure];

    thisWfWidthOSIPlot = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);
    subplot(1,3,1)
    scatter(thiswfWidth,thisTable.pure,10,'k','filled')
    xlabel('Waveform width (ms)')
    ylabel('Pure OSI')
    ylim([0 1])
    legend off;
    ax = gca;
    xticks(0:0.2:1.6)
    yticks(0:0.2:1)
    ax.FontSize = 14;
    
    subplot(1,3,2)
    scatter(thiswfWidth,thisTable.noOpto,10,'k','filled')
    xlabel('Waveform width (ms)')
    ylabel('Interleaved noOpto OSI')
    ylim([0 1])
    legend off;
    ax = gca;
    xticks(0:0.2:1.6)
    yticks(0:0.2:1)
    ax.FontSize = 14; 
    
    subplot(1,3,3)
    scatter(thiswfWidth,thisTable.opto,10,'k','filled')
    xlabel('Waveform width (ms)')
    ylabel('Interleaved opto OSI')
    ylim([0 1])
    legend off;
    ax = gca;
    xticks(0:0.2:1.6)
    yticks(0:0.2:1)
    ax.FontSize = 14;

    sgtitle(fieldName)
    fPath = fullfile(cd,'combined_analysis');
    try
        saveFigurePdf(thisWfWidthOSIPlot, fullfile(fPath, append(fieldName,'_wfWidthOSI')))    
    catch ME
        disp(['failed to save figure'])
    end
end

combinedWfWidthOSIPlots = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);
subplot(2,2,1)
scatter(combinedwfWidth,combinedPure,10,'k','filled')
xlabel('Waveform width (ms)')
ylabel('Pure OSI')
ylim([0 1])
legend off;
ax = gca;
xticks(0:0.2:1.6)
yticks(0:0.2:1)
set(gca,'LineWidth',2,'FontSize',14)

subplot(2,2,2)
scatter(combinedwfWidth,combinedNoOpto,10,'k','filled')
xlabel('Waveform width (ms)')
ylabel('Interleaved noOpto OSI')
ylim([0 1])
legend off;
ax = gca;
xticks(0:0.2:1.6)
yticks(0:0.2:1)
set(gca,'LineWidth',2,'FontSize',14)

subplot(2,2,3)
scatter(combinedwfWidth,combinedOpto,10,'k','filled')
xlabel('Waveform width (ms)')
ylabel('Interleaved opto OSI')
ylim([0 1])
legend off;
ax = gca;
xticks(0:0.2:1.6)
yticks(0:0.2:1)
set(gca,'LineWidth',2,'FontSize',14)
box off;


% calculate OSI average for units < 0.45 ms and units > 0.45
FS_idx = find(combinedwfWidth < 0.45);
RS_idx = find(combinedwfWidth >= 0.5);

pure_FS = combinedPure(FS_idx);
pure_RS = combinedPure(RS_idx);
opto_FS = combinedOpto(FS_idx);
opto_RS = combinedOpto(RS_idx);

categories = {'FS pure','RS pure','FS opto','RS opto'};

subplot(2,2,4)
swarmchart(repmat(1,1,length(pure_FS)),pure_FS,10,'k','filled');
hold on;
scatter(1,nanmean(pure_FS),40,'r','filled')
hold on;

swarmchart(repmat(2,1,length(pure_RS)),pure_RS,10,'k','filled');
hold on;
scatter(2,nanmean(pure_RS),40,'r','filled')
hold on;

swarmchart(repmat(3,1,length(opto_FS)),opto_FS,10,'k','filled');
hold on;
scatter(3,nanmean(opto_FS),40,'r','filled')
hold on;

swarmchart(repmat(4,1,length(opto_RS)),opto_RS,10,'k','filled');
xticks(1:4)
scatter(4,nanmean(opto_RS),40,'r','filled')
hold on;

xticklabels({'Pure FS','Pure RS','Opto FS','Opto RS'})
ax = gca;
set(gca,'LineWidth',2,'FontSize',14)
box off;

pVal_pure = ranksum(pure_FS,pure_RS);
pVal_opto = ranksum(opto_FS,opto_RS);
pVal_FS = signrank(pure_FS,opto_FS);
pVal_RS = signrank(opto_RS,opto_RS);

pValString = append('p_pure = ',num2str(pVal_pure),...
    ' p_opto = ',num2str(pVal_opto),...
    ' p_FS = ',num2str(pVal_FS),...
    ' p_RS = ',num2str(pVal_RS));

annotation('textbox', [0.9, 0.9, 0.1, 0.05], 'String', pValString, ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center');

sgtitle('Combined waveform width')
saveFigurePdf(combinedWfWidthOSIPlots, fullfile(fPath, append(fieldName,'_combinedWfWidthOSI')))    
close all;


%% OMI sanity
fields = fieldnames(OMI_struct);

for i = 1:length(fields)
    fieldName = fields{i};
    thisOMITable = OMI_struct.(fieldName).optoModulation_Table;
    
    thisOMIPlot = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
    scatter(thisOMITable.OMI,thisOMITable.otherOMI,10,'k','filled')
    hold on;
    plot(-1:0.1:1,-1:0.1:1,'r-')
    xlabel('Opto modulation index')
    ylabel('Alternative opto modulation index')
    xlim([-1 1])
    ylim([-1 1])
    legend off;
    ax = gca;
    xticks(-1:0.2:1)
    yticks(-1:0.2:1)
    ax.FontSize = 14;

    sgtitle(fieldName)
    fPath = fullfile(cd,'combined_analysis');
    try
        saveFigurePdf(thisOMIPlot, fullfile(fPath, append(fieldName,'_OMISanityCheck')))    
    catch ME
        disp(['failed to save figure'])
    end
end 

close all;



%% OMI vs waveform width
fields = fieldnames(wfWidth_struct);

combinedwfWidth = [];
combinedOMI = [];

for i = 1:length(fields)
    fieldName = fields{i};
    thisTable = OMI_struct.(fieldName).optoModulation_Table;
    thiswfWidth = wfWidth_struct.(fieldName).waveformWidth(2,:);
    combinedwfWidth = [combinedwfWidth thiswfWidth];
    combinedOMI = [combinedOMI; thisTable.OMI];

    thisWfWidthOMIPlot = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
    scatter(thiswfWidth,thisTable.OMI,10,'k','filled')
    xlabel('Waveform width (ms)')
    ylabel('Opto-modulation index')
    ylim([-1 1])
    legend off;
    ax = gca;
    xticks(0:0.2:1.6)
    yticks(-1:0.2:1)
    ax.FontSize = 14;

    sgtitle(fieldName)
    fPath = fullfile(cd,'combined_analysis');
    try
        saveFigurePdf(thisWfWidthOMIPlot, fullfile(fPath, append(fieldName,'_wfWidthOMI')))    
    catch ME
        disp(['failed to save figure'])
    end
end

combinedWfWidthOMIPlot = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
scatter(combinedwfWidth,combinedOMI,10,'k','filled')
xlabel('Waveform width (ms)')
ylabel('Opto-modulation index')
ylim([-1 1])
legend off;
ax = gca;
xticks(0:0.2:1.6)
yticks(-1:0.2:1)
ax.FontSize = 14;
close all;

%% OMI vs depth
fields = fieldnames(OMI_struct);

combinedOMI = [];
combinedwfDepth = [];
combinedwfWidth = [];
for i = 1:length(fields)
    fieldName = fields{i};
    thisOMITable = OMI_struct.(fieldName).optoModulation_Table;
    thiswfDepthTable = wfDepth_struct.(fieldName).waveformDepth_Table;
    thiswfWidth = wfWidth_struct.(fieldName).waveformWidth(2,:);
    combinedwfWidth = [combinedwfWidth thiswfWidth];
    combinedOMI = [combinedOMI; thisOMITable.OMI];
    combinedwfDepth = [combinedwfDepth; thiswfDepthTable.yCoord];

    thisOMIDepthPlot = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
    scatter(thisOMITable.OMI,thiswfDepthTable.yCoord,10,'k','filled');
    xlabel('Opto modulation index')
    ylabel('Depth from deepest part of probe (um)')
    xlim([-1 1])
    legend off;
    ax = gca;
    xticks(-1:0.2:1)
    ax.FontSize = 14;

    sgtitle(fieldName)
    fPath = fullfile(cd,'combined_analysis');
    try
        saveFigurePdf(thisOMIDepthPlot, fullfile(fPath, append(fieldName,'_OMIDepth')))    
    catch ME
        disp(['failed to save figure'])
    end
end 


% only plot FS
FS_idx = find(combinedwfWidth < 0.45); 
RS_idx = find(combinedwfWidth > 0.45);

scatter(combinedOMI(FS_idx),combinedwfDepth(FS_idx),10,'k','filled')



scatter(combinedOMI(RS_idx),combinedwfDepth(RS_idx),10,'k','filled')


combinedOMIPlot = figure('Renderer', 'painters', 'Position', [10 10 800 800]);
histogram(combinedOMI,-1:0.05:1,'FaceColor','black');
hold on;
xline(mean(combinedOMI),'r-','LineWidth',2.0)
xlabel('Opto-modulation index')
ylabel('Counts')
box off;
ax = gca;
ax.FontSize = 14;
sgtitle('Combined OMI')
saveFigurePdf(combinedOMIPlot, fullfile(fPath, append(fieldName,'_combinedOMI')))    
close all;