function [amplitude,pdf,fraction_missing] = percSpikesMissing(gsaparams)
%
% Written by Chris (Seong Yeol) An
% san11@jh.edu
% 8/21/2024
%
% https://spikeinterface.readthedocs.io/en/stable/modules/qualitymetrics/amplitude_cutoff.html
%
% ------
% Inputs
% ------
% gsaparams: parameter structure to run getSpikeAmplitudes
%
% ------
% Outputs
% ------
% amplitude: array of size 1 x spikeCounts
% pdf: pdf of gaussian fit to the amplitude histogram
% fraction_missing: 0.5 if more than 0.5, but otherwise represents a
% fraction of spikes that are missing from the gaussian assumption.


histogram_smoothing_value = 3;
spikeAmplitudeThisUnit = getSpikeAmplitudes(gsaparams);
numOfBins = 100;

amplitude = (spikeAmplitudeThisUnit.amplitude * 0.5 / 8192) / 80 * 1000000;

h = histfit(spikeAmplitudeThisUnit.amplitude,numOfBins,'normal')

[N,edges] = histcounts(spikeAmplitudeThisUnit.amplitude,numOfBins,'Normalization','probability');
support = edges(1:end-1) + diff(edges)/2;

% Gaussian smoothing?
pdf = imgaussfilt(N,histogram_smoothing_value);
[~, peak_index] = max(pdf); % largest value of the smoothed gaussian

% Find the index where PDF becomes closest to initial value
[~,G] = min(abs(pdf(peak_index:end) - pdf(1)));
G = G + peak_index - 1;

% calculate bin size
bin_size = mean(diff(support));

% estimate the fraction of missing spikes
fraction_missing = sum(pdf(G:end)) * bin_size;

fraction_missing = min(fraction_missing,0.5)

end