# MEA_DOI

## Project Overview
This repository contains MATLAB scripts and tools for analyzing multi-electrode array (MEA) data with a focus on studying the effects of DOI (2,5-dimethoxy-4-iodoamphetamine), a psychedelic compound, on neural activity. This project is part of the Psychoplastogens research at Hai's Lab at UW Madison.

## Structure
- **TDTSDK/**: TDT (Tucker-Davis Technologies) MATLAB SDK for analyzing neurophysiological data
  - **APIStreamer/**: Tools for streaming data from TDT hardware
  - **OpenExLive/**: Live data acquisition tools
  - **SynapseAPI/**: API for interacting with TDT Synapse software
  - **SynapseLive/**: Live data analysis tools for Synapse
  - **TDTbin2mat/**: Converts TDT binary files to MATLAB format
  - **UDP/**: Network communication tools

- **scripts/**: Analysis scripts for MEA data processing
  - **Psychoplastogen_FFT_1D.m**: Frequency analysis of MEA recordings
  - **histo_Rasters_dataset.m**: Generates histograms and raster plots of neural activity
  - **raster_auto.m**: Automated raster plot generation
  - **spectro.m** and **spectrogram_dev.m**: Spectral analysis tools

- **docs/**: Documentation and reference papers
  - **OfflineAnalysisToolsMatlab.pdf**: Guide for offline data analysis
  - **Sabanovic_2022_Structural_and_behavioural.pdf**: Reference research paper

## Key Features
- Preprocessing of neural recordings with bandpass and notch filters
- Spike detection and analysis of neural activity
- Raster plot and histogram generation for visualizing activity
- Comparison of baseline vs. drug-exposed neural recordings
- Spectral analysis (FFT) of neural signals

## Usage
The scripts process TDT Synapse 'tank' recordings, applying appropriate filters to detect and analyze neural spikes. Data analysis focuses on comparing control/baseline recordings with recordings after DOI application, looking at metrics such as:
- Firing frequency
- Number of active channels
- Total spike counts
- Spectral characteristics

Results are saved as both data (.csv) and visualizations (.png).

## Dependencies
- MATLAB
- TDT MATLAB SDK (included in repository)
- TDT Synapse Software (for data acquisition)