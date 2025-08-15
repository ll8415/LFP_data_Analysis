# LFP data Analysis GUI

## Overview
This project provides a graphical user interface (GUI) for loading, visualizing, and analyzing local field potential (LFP) signals.  
It supports waveform plotting, power spectral density (PSD) analysis, continuous wavelet transforms (CWT), and multi-band signal visualization.

## Features
- Load LFP data in Neo-compatible formats  
- Time-domain signal visualization  
- Power Spectral Density (PSD) computation  
- Continuous Wavelet Transform (CWT) analysis  
- Multi-band filtering and visualization  
- Channel selection and time range control  

## Requirements
- Python 3.10+  
- Dependencies:
  ```bash
  pip install neo numpy scipy matplotlib quantities pywt
  ```
<img width="1778" height="1218" alt="2025-08-15 164345" src="https://github.com/user-attachments/assets/94c1e36d-151f-4fd4-afe4-d334680791dd" />



## Notes
- Ensure your LFP file is in a Neo-supported format (e.g., `.ns5`, `.nix`, `.h5`).  
- `t_start` and `t_stop` must be within the duration of the recorded signal.  
- Use the console output to check for errors when loading or plotting data.  
