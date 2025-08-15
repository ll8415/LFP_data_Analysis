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

## How to Run
1. Save the script (`lfp_gui_runner_compat_v6_en.py`) in your working directory.  
2. Open a terminal and run:
   ```bash
   python lfp_gui_runner_compat_v6_en.py
   ```
3. The GUI window will appear.

## Usage
1. **Load Data**: Click "Open" and select your LFP file.  
2. **Select Channel**: Choose the desired channel from the dropdown.  
3. **Plot Signal**: Click "Plot Signal" to visualize the LFP waveform.  
4. **PSD Analysis**: Click "Plot PSD" to see the frequency-domain representation.  
5. **Wavelet Analysis**: Select a time point and click "Plot CWT".  
6. **Multi-band Analysis**: Choose start time, duration, and channel; click "Plot Multi-band".  

## Notes
- Ensure your LFP file is in a Neo-supported format (e.g., `.ns5`, `.nix`, `.h5`).  
- `t_start` and `t_stop` must be within the duration of the recorded signal.  
- Use the console output to check for errors when loading or plotting data.  
