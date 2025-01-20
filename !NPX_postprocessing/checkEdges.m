%% Check edge counts
% Here, I want to look at CatGT edge extracted files and compare the number of
% entries with the RPi file. 
% This makes sure that CatGT edge is accurately extracting the information.
 
% NOTES ON 8/6/2024
% COMPARING LICKS AND WATER DELIVERY IS VERY DIFFICULT DUE TO HOW THE DATA
% WAS ACQUIRED (AND SYNCHRONIZED) BETWEEN NIDAQ AND THE RPI. THERE SHOULD
% BE A MECHANISM ON SOMEHOW TELLING THE RPI TO START RECORDING LICKS,
% ETC... MAYBE ADD A PIECE OF LINE ONTO THE CODE IN RPI SO THAT WHEN YOU
% PRESS X OR SOMETHING IT PRINTS A LINE "START ACQUISITION" AND "END
% ACQUISITION..."

%% Set up paths to import the data
% Extract session names from RPi files

parentFolder = 'D:\NPX_data\DATA'; % Parent folder for holding all NPX data
RPiFolder = fullfile(parentFolder,'RPi_logs'); % Folder for holding RPi files
RPiLogList = dir(RPiFolder); % list of file names
RPiSessionNames = {}; % This will hold the sessionNames from the RPi log files
counter = 0;
for i = 1:length({RPiLogList.name})
    temp = RPiLogList(i).name;
    tempSplit = strsplit(temp,'.');
    if ~isempty(tempSplit{1,1})
        counter = counter + 1;
        RPiSessionNames{counter} = tempSplit{1,1};
    end
end

% Find matching session names using the RPi files
recordingDirs = dir(parentFolder);
dirFlags = [recordingDirs.isdir];
recordingFolders = {recordingDirs(dirFlags).name}; % Only extract folders
sessionFlags = contains(recordingFolders,RPiSessionNames); % Find sessions that are in RPi folder
recordingSessionNames = recordingFolders(sessionFlags); % This will hold the matching NPX


%% Now we want to iterate through groups of recordingSessionNames and group the CatGT edge files together

for i = 2:length(RPiSessionNames)
    thisRPiSessionName = RPiSessionNames{i};
    thisRecordingSessionNames = recordingSessionNames(contains(recordingSessionNames,thisRPiSessionName));
    
    combinedCatGTVisStim = []; % This will hold all instances of visual stim
    combinedCatGTLicks = []; % This will hold all instances of licks
    combinedCatGTWater = []; % This will hold all instances of water
    % Now we go into subfolders of thisRecordingSessionNames
    for j = 1:length(thisRecordingSessionNames)
        thisRecordingSubFolder = fullfile(parentFolder,thisRecordingSessionNames{j},thisRecordingSessionNames{j});

        % Load VisStim
        thisCatGTVisStimFile = fullfile(thisRecordingSubFolder,'*_1500.txt');
        thisCatGTVisStimDataFile = dir(thisCatGTVisStimFile);
        thisCatGTVisStim = readmatrix(fullfile(thisCatGTVisStimDataFile.folder,thisCatGTVisStimDataFile.name)); % This loads the CatGT file
        combinedCatGTVisStim = [combinedCatGTVisStim thisCatGTVisStim'];

        % Find where one set of vis stim ends
        idx = find(diff(thisCatGTVisStim) >= 4 == 1);

        % Load licks
        thisCatGTLicksFile = fullfile(thisRecordingSubFolder,'*_4_0.txt');
        thisCatGTLicksDataFile = dir(thisCatGTLicksFile);
        thisCatGTLicks = readmatrix(fullfile(thisCatGTLicksDataFile.folder,thisCatGTLicksDataFile.name)); % This loads the CatGT file
        thisCatGTLicks_truncated = thisCatGTLicks(thisCatGTLicks >= thisCatGTVisStim(idx) & thisCatGTLicks <= thisCatGTVisStim(idx+1)); % Remove any element that happens before first vis stim.

        combinedCatGTLicks = [combinedCatGTLicks thisCatGTLicks_truncated'];

        % Load water
        thisCatGTWaterFile = fullfile(thisRecordingSubFolder,'*_7_0.txt');
        thisCatGTWaterDataFile = dir(thisCatGTWaterFile);
        thisCatGTWater = readmatrix(fullfile(thisCatGTWaterDataFile.folder,thisCatGTWaterDataFile.name)); % This loads the CatGT file
        thisCatGTWater_truncated = thisCatGTWater(thisCatGTWater >= thisCatGTVisStim(1) & thisCatGTWater <= thisCatGTVisStim(end)); % Remove any element that happens before first vis stim.

        combinedCatGTWater = [combinedCatGTWater thisCatGTWater_truncated'];
    end
        
    % Now we want to compare the combined files and RPi 
    loadRPiSession = fullfile(RPiFolder,append(thisRPiSessionName,'.txt'));
    RPiSessionData = readtable(loadRPiSession);
    % Row where visual stimulation first appears
    firstVisStimIdx = find(diff(contains(RPiSessionData.Var1,'map'))==-1,1);
    lastVisStimIdx = find(diff(contains(RPiSessionData.Var1(firstVisStimIdx:end),'map'))==1,1) + firstVisStimIdx;
    RPiSessionData_truncated = RPiSessionData(firstVisStimIdx:lastVisStimIdx,:);
    counts_RPiLick = sum(matches(RPiSessionData_truncated.Var1,'Lick'));
    counts_RPiWater = sum(matches(RPiSessionData_truncated.Var1,'MANUAL WATER DELIVERED'));
    counts_RPiVisSineMap = sum(matches(RPiSessionData_truncated.Var1,'Sine map'));
    counts_RPiVisGratingMap = sum(matches(RPiSessionData_truncated.Var1,'Grating map'));
    counts_RPiVisSquareMap = sum(matches(RPiSessionData_truncated.Var1,'Square map'));
    counts_RPiVisMap = counts_RPiVisSquareMap + counts_RPiVisSineMap + counts_RPiVisGratingMap;

    if counts_RPiVisMap == length(combinedCatGTVisStim)
        message = append(num2str(length(combinedCatGTLicks)),'vs.',num2str(counts_RPiLick),'   ',...
            num2str(length(combinedCatGTWater)),'vs.',num2str(counts_RPiWater));
        f = msgbox(message);
        uiwait(msgbox(append(thisRPiSessionName,' Length of RPi visual stim stamps and CatGT stamps are the same!')));
    else
        f = msgbox(message)
        uiwait(msgbox(append(thisRPiSessionName,' Length of RPi visual stim stamps and CatGT stamps are NOT the same!')));
    end
end




counts_RPiWater = sum(matches(RPiSessionData_truncated.Var1,'MANUAL WATER DELIVERED'));
