%% Piece of code to combine PDFS generated from functions in master_analysis

% This should be the parent folder that has all subfolders with PDFs
rootFolder = cd; % usually kilosort_output folder

% Search for all pdfs in the subfolders
pdfFiles = dir(fullfile(rootFolder,'**','*.pdf'));

% Initialize a map to hold filenames and their corresponding file paths
pdfMap = containers.Map();

% Populate the map
for k = 1:length(pdfFiles)
    [~,fileName,ext] = fileparts(pdfFiles(k).name);
    fullPath = fullfile(pdfFiles(k).folder, pdfFiles(k).name);
    key = strcat(fileName,ext);

    if isKey(pdfMap,key)
         pdfMap(key) = [pdfMap(key), {fullPath}];
    else
        pdfMap(key) = {fullPath};
    end
end


% Specify the output folder for combined PDFs
outputFolder = 'combined_pdf';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Combine PDFs with the same names
keys = pdfMap.keys;
for i = 1:length(keys)
    key = keys{i};
    pdfList = pdfMap(key);

    % Custom sort the order of PDFs
    customOrder = {'waveforms','Orientation','OptoOrientation','PolarOrientation','OptoScatter'};
    folderPaths = cellfun(@(x) fileparts(x), pdfList, 'UniformOutput', false);
    [~,folderNames,~] = cellfun(@(x) fileparts(x), folderPaths, 'UniformOutput', false);

    %[~,sortIdx] = ismember(folderNames, customOrder);
    sortedPdfList = cell(size(customOrder));
    for i = 1:length(customOrder)
        % Find the index in folderNames that matches the current customOrder
        idx = find(strcmp(folderNames, customOrder{i}));
        
        if ~isempty(idx)
            sortedPdfList{i} = pdfList{idx};
        else
            warning('Folder %s not found in pdfList.', customOrder{i});
        end
    end    


    % Merge the sorted PDFs
    outputFilePath = fullfile(outputFolder, key);
    
    % Use the fex-pdfjoin or a similar function to join the PDFs
    % Example using pdfjoin (install from File Exchange)
    mergePdfs(sortedPdfList, outputFilePath);
end