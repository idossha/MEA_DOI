% 1. Housekeeping
close all; clear; clc;

% 2. Set initial parameters
SDKPATH = '/Volumes/ahai/TDT/TDTMatlabSDK';

% DOI path:
BASEPATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/NewBatch_09140223';
desktopPath = '/Users/idohaber/Desktop/Final_DOI'; % Set this to your output directory
addpath(genpath(SDKPATH));

% Parameters for data processing
STORE = 'Wav1';
CHANNELS = 1:64; % channels we want to process for raster
BP_FILTER = [300 2500]; % Band-pass filter parameters
f0 = 60; % Notch filter center frequency
fs = 48000; % Sampling frequency
fn = fs/2;
freqRatio = f0/fn;
[bw, a] = iirnotch(freqRatio, freqRatio/35);

% For quick testing
datasets = {'IdoControl-230914-130200_#1'};

% Initialize firing frequency storage
firingFrequencyPerDataset = zeros(length(datasets), 1);

% Loop through each dataset
for i = 1:length(datasets)
    ENDPATH = datasets{i};
    BLOCKPATH = fullfile(BASEPATH, ENDPATH);

    datasetSpikes = 0;
    datasetActiveChannels = 0;
    recordLengthSecs = 0;
    spikeCountsPerChannel = zeros(length(CHANNELS), 1);

    rasterData = cell(length(CHANNELS), 1);
    for ch = CHANNELS
        data = TDTbin2mat(BLOCKPATH, 'STORE', STORE, 'CHANNEL', ch);
        dataFiltered = TDTdigitalfilter(data, STORE, BP_FILTER, 'ORDER', 5);
        dataFiltered.streams.(STORE).data = filtfilt(bw, a, dataFiltered.streams.(STORE).data);
        
        spikes = TDTthresh(dataFiltered, STORE, 'MODE', 'auto', 'POLARITY', -1, 'STD', 6.5, 'TAU', 5);
        rasterData{ch} = spikes.snips.Snip.ts;

        if ch == CHANNELS(i)
            fs = data.streams.(STORE).fs;
            nSamples = length(data.streams.(STORE).data);
            recordLengthSecs = nSamples / fs;
        end

        if ~isempty(rasterData{ch})
            spikeCount = length(rasterData{ch});
            datasetSpikes = datasetSpikes + spikeCount;
            spikeCountsPerChannel(ch) = spikeCount;
            datasetActiveChannels = datasetActiveChannels + 1;
        end
    end
    
    % Calculate and store firing frequency for the current dataset
    firingFrequencyPerDataset(i) = datasetSpikes / recordLengthSecs;
   
    % Plotting Raster Plot
    figure; hold on;
    for ch = 1:length(CHANNELS)
        spikeTimes = rasterData{ch};
        for spike = 1:length(spikeTimes)
            plot([spikeTimes(spike), spikeTimes(spike)], [ch-0.4, ch+0.4], 'k');
        end
    end
    title(sprintf('%s\nActive Channels = %d, Total Spikes = %d, Firing Frequency = %.2f spikes/s', ENDPATH, datasetActiveChannels, datasetSpikes, firingFrequencyPerDataset(i)));
    xlabel('Time (s)'); ylabel('Channel'); axis tight;
    saveas(gcf, fullfile(desktopPath, sprintf('raster_%s.png', ENDPATH)), 'png');
    close;

    % Spike Count per Channel Plot
    figure;
    bar(1:length(CHANNELS), spikeCountsPerChannel, 'FaceColor', 'b');
    title(sprintf('Spike Count per Channel for %s', ENDPATH));
    xlabel('Channel Number'); ylabel('Spike Count');
    grid on; axis tight;
    saveas(gcf, fullfile(desktopPath, sprintf('spike_count_per_channel_%s.png', ENDPATH)), 'png');
    close;

end


%%

TDTfft(data.streams.Wav1, 2, 'FREQ', [0, 100], 'NUMAVG', 20,'SPECPLOT', 1);