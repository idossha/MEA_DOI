% 1. Housekeeping
close all; clear; clc;
% 2. Set initial parameters
SDKPATH = '/Volumes/ahai/TDT/TDTMatlabSDK';

% New Batch path:
%BASEPATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/NewBatch_09140223';

%Ketanserin path:

BASEPATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/Ketanserin_11012023';

% datasets = {...
%     'IdoControl-230914-130200_#1', ...
%     'IdoControl-230914-131548_#2', ...
%     'IdoControl-230914-132855_#3', ...
%     'IdoControl-230914-154022_#4', ...
%     'IdoControl-230914-155318_#5', ...
%     'IdoControl-230914-160601_#6', ...
%     'IdoDOI-230914-142502_#1', ...
%     'IdoDOI-230914-143740_#2', ...
%     'IdoDOI-230914-144945_#3', ...
%     'IdoDOI-230914-161838_#4', ...
%     'IdoDOI-230914-163140_#5', ...
%     'IdoDOI-230914-164512_#6'};

% control1:
% datasets = {'IdoControl-230914-130200_#1','IdoControl-230914-131548_#2','IdoControl-230914-132855_#3','IdoControl-230914-154022_#4','IdoControl-230914-155318_#5','IdoControl-230914-160601_#6'};
% DOI1:
% datasets = {'IdoDOI-230914-144945_#3','IdoDOI-230914-161838_#4','IdoDOI-230914-163140_#5','IdoDOI-230914-164512_#6'};

% control2:
%datasets = {'IdoControl-231101-144046_#1','IdoControl-231101-145330_#2','IdoControl-231101-150512_#3'};
% Ketanserin
 datasets = {'IdoKetanserin-231101-133725_#1','IdoKetanserin-231101-134946_#2','IdoKetanserin-231101-140237_#3'};

addpath(genpath(SDKPATH));

% Parameters for data processing
STORE = 'Wav1';
CHANNELS = 1:64; % channels we want to process for raster
BP_FILTER = [300 2500]; % Band-pass filter parameters
f0 = 60; % Notch filter center frequency
fs = 48000; % Sampling frequency
fn = fs/2;
freqRatio = f0/fn;
[bw, a] = iirnotch(freqRatio, freqRatio/0.5);

desktopPath = '/Users/idohaber/Desktop/PP_Output'; % Set this to your output directory

for i = 1:length(datasets)
    ENDPATH = datasets{i};
    BLOCKPATH = fullfile(BASEPATH, ENDPATH);

    % Initialize the counters for each dataset
    datasetSpikes = 0; % Total spikes in the current dataset
    datasetActiveChannels = 0; % Active channels in the current dataset
    recordLengthSecs = 0; % Initialize recording length in seconds for the current dataset

    % 5. Data processing and figure generation
    rasterData = cell(length(CHANNELS), 1);
    for ch = CHANNELS
        data = TDTbin2mat(BLOCKPATH, 'STORE', STORE, 'CHANNEL', ch);
        dataFiltered = TDTdigitalfilter(data, STORE, BP_FILTER, 'ORDER', 5);
        dataFiltered.streams.(STORE).data = filtfilt(bw, a, dataFiltered.streams.(STORE).data);
        
        spikes = TDTthresh(dataFiltered, STORE, 'MODE', 'auto', 'POLARITY', -1, 'STD', 6, 'TAU', 5);
        rasterData{ch} = spikes.snips.Snip.ts;

        % For recording length, ensure it's calculated once per dataset
        if ch == CHANNELS(1)
            fs = data.streams.(STORE).fs;
            nSamples = length(data.streams.(STORE).data);
            recordLengthSecs = nSamples / fs;
        end

        % Update spike counts and active channel counts
        if ~isempty(rasterData{ch})
            datasetSpikes = datasetSpikes + length(rasterData{ch});
            datasetActiveChannels = datasetActiveChannels + 1;
        end
    end

    % 6. Plotting Raster Plot
    figure; hold on;
    for ch = 1:length(CHANNELS)
        spikeTimes = rasterData{ch};
        for spike = 1:length(spikeTimes)
            plot([spikeTimes(spike), spikeTimes(spike)], [ch-0.4, ch+0.4], 'k');
        end
    end
    
    % Calculate overall firing frequency for the current dataset
    firingFrequency = datasetSpikes / recordLengthSecs;

    % Printing the results for the current dataset
    fprintf('Dataset %s: Active Channels = %d, Total Spikes = %d\n', ENDPATH, datasetActiveChannels, datasetSpikes);
    fprintf('Overall Firing Frequency: %f spikes per second\n', firingFrequency);
    xlabel('Time (s)'); ylabel('Channel'); title(sprintf('%s\nActive Channels = %d, Total Spikes = %d, Firing Frequency = %.2f spikes/s', ENDPATH, datasetActiveChannels, datasetSpikes, firingFrequency)); axis tight;
    
    % Save the figure
    filename = sprintf('raster_%s.png', ENDPATH);
    fullFilePath = fullfile(desktopPath, filename);
    saveas(gcf, fullFilePath, 'png');
    close; % Close the figure to save memory
end
