function presenceRatio = calculate_presence_ratio(spikeTimes,binSize)
% presenceRatio : fraction of bins (of bin size presenceRatioBinSize) that
% contain at least 5% of the largest spike count per bin

    % divide recordings times in chunks
    presenceRatio_bins = min(spikeTimes):binSize:max(spikeTimes);
    % count number of spikes in each chunk 
    spikesPerBin = histcounts(spikeTimes,presenceRatio_bins);
    fullBins = spikesPerBin >= 0.05*prctile(spikesPerBin, 90); % calculate the 90th percentile
    presenceRatio = sum(fullBins)/length(spikesPerBin);

end