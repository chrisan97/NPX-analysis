function dataStruct = searchAndLoad(currentDir, targetFileName, dataStruct, rootDir)
    % Get list of all files and folders in the current directory
    contents = dir(currentDir);

    % Loop through each item in the directory
    for i = 1:length(contents)
        item = contents(i);
        
        % Skip the '.' and '..' folders
        if strcmp(item.name, '.') || strcmp(item.name, '..')
            continue;
        end
        
        % Construct full path to the item
        fullPath = fullfile(currentDir, item.name);
        
        % Check if the item is a folder
        if item.isdir
            % Check if the folder is named "analysis_variables"
            if strcmp(item.name, 'analysis_variables')
                % Look for the specified .mat file in this folder
                matFile = dir(fullfile(fullPath, [targetFileName, '.mat']));
                
                if ~isempty(matFile)
                    % Load the .mat file and store it in the structure
                    matFilePath = fullfile(fullPath, matFile.name);
                    disp(['Loading ', matFilePath]);
                    try
                        loadedData = load(matFilePath);
                        % Use the names of the immediate subfolders for unique keys
                        parentFolder = fileparts(fullPath);
                        parentFolderName = strrep(parentFolder, rootDir, '');
                        splitFolderName = strsplit(parentFolderName,"\");
                        parentFolderName = splitFolderName{2};
                        splitFolderName = strsplit(parentFolderName,'-');
                        parentFolderName = [splitFolderName{2} '_' splitFolderName{1}];
                        subfolderKey = strcat(parentFolderName);
                        dataStruct.(subfolderKey) = loadedData;
                    catch ME
                        warning('Failed to load %s: %s', matFilePath, ME.message);
                    end
                    return; % Stop searching after loading the file
                else
                    disp(['No .mat files named ', targetFileName, ' found in ', fullPath]);
                end
            else
                % Recursively search in the subfolder
                dataStruct = searchAndLoad(fullPath, targetFileName, dataStruct, rootDir);
            end
        end
    end
end