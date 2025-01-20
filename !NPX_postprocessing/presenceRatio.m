function presence_ratio = presenceRatio(spikeTimes,binSize)
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/21/2024
%
% https://spikeinterface.readthedocs.io/en/stable/modules/qualitymetrics/presence_ratio.html
%
% ------
% Inputs
% ------
% spikeTimes: spike timestamps in seconds
% binSize: binSize to use to assess presence_ratio (this should also be in
% seconds)
%
% ------
% Outputs
% ------
% presence_ratio: fraction of bins that contain at least 5% of the largest
% spike count per bin
%

presence_ratio_bins = spikeTimes(1):binSize:spikeTimes(end); 
spikesPerBin = arrayfun(@(x) sum(spikeTimes >= presence_ratio_bins(x) & spikeTimes < presence_ratio_bins(x+1)),1:length(presence_ratio_bins)-1);
fullbins = spikesPerBin >= 0.05 * prctile(spikesPerBin,90);
presence_ratio = sum(fullbins)/length(spikesPerBin);


end