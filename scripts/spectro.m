% 1. Housekeeping
close all; clear; clc;

% 2. Set initial parameters
SDKPATH = '/Volumes/ahai/TDT/TDTMatlabSDK';
PATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/NewBatch_09140223/IdoControl-230914-130200_#1';

% Load data using TDTbin2mat
data = TDTbin2mat(PATH, 'T1', 0 , 'T2', 500 );

% Assuming you have continuous data in a stream called 'Wav1'
fs = data.streams.Wav1.fs; % Sampling frequency
signal = data.streams.Wav1.data; % Signal data

% Bandpass filter parameters
bpFilt = designfilt('bandpassiir', 'FilterOrder', 20, 'HalfPowerFrequency1', 300, 'HalfPowerFrequency2', 3000, 'SampleRate', fs);

% Apply bandpass filter
filteredSignal = filtfilt(bpFilt, signal);

% Notch filter to remove 60 Hz line noise
d = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', fs);

% Apply notch filter
filteredSignal = filtfilt(d, filteredSignal);

t = (0:length(signal)-1)/fs; % Time vector
%%
% Plot raw and filtered signal
figure;
subplot(2,1,1);
plot(t, signal(60,:));
title('Raw Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, filteredSignal(60,:));
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');
%%

% Create a structure for filtered data mimicking 'Wav1' structure
filteredData.streams.Wav1.data = filteredSignal;
filteredData.streams.Wav1.fs = fs;

% Call TDTfft on the filtered data instead of the raw data
TDTfft(filteredData.streams.Wav1, 5, 'NUMAVG', 20, 'SPECPLOT', 1, 'FREQ', [400, 2000]);

set(gcf, 'Position', [100 100 750 800])
