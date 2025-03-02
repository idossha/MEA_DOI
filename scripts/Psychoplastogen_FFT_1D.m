%= ReadTankMEA.m | Ilhan Bok | Hai Lab 2021
%= MATLAB Script
%= Last updated: 16 Mar. 2021

%= Purpose: Reads TDT Synapse 'tank' recordings and displays onto a
%=          labeled plot

% Clean up previous runs
clear
clc
close all

% 2. Set initial parameters
SDKPATH = '/Volumes/ahai/TDT/TDTMatlabSDK';


%DOI path:
BASEPATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/NewBatch_09140223/IdoControl-230914-130200_#1';

desktopPath = '/Users/idohaber/Desktop/Final_DOI'; % Set this to your output directory
addpath(genpath(SDKPATH));


%= Everything below this comment is automatic. Contact Ilhan for
%= assistance.

%STORE = 'LFP1';
%CHANNEL = 1;
%T2 = 10;

disp('Performing Bandpass of [300 3051] Hz...');

% Read out data from Tank
data = TDTbin2mat(BASEPATH, 'TYPE', {'epocs', 'scalars', 'streams'});
data = TDTdigitalfilter(data, 'Wav1', [300 3051], 'ORDER', 1); % Max 3051.75 per Nyquist
%all_streams = highpass(data.streams.Wav1.data',0.0491)';
all_streams= data.streams.Wav1.data;
%%
disp('Plotting...');
channel=(1:64);
for i = (1)
    Y = fft(all_streams(channel(i),:));
    
    Fs = 200e3;              % Sampling frequency                    
    T = 1/Fs;               % Sampling period       
    L = length(all_streams(channel(i),:));  % Length of signal
    
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f( 1:find( f > 1000, 1 ) ),P1( 1:find( f > 1000, 1 ) ))
    hold on
end
ylim([0 2.5e-7])
set(0, 'defaultFigureRenderer', 'painters');
set(gcf, 'Position', get(0, 'Screensize'));
%set(gcf,'position',[0 0 595 776])
saveas(gcf,strcat(files_to_read{l},'_1000Hz_FFT.png'),'png');
close all
