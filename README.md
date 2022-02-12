# Posdoc
## Description
Files, data and programs useful for posdoctoral project and paper.

## Content:
1. journal.org
  - Contains tables, to do suff and other useful material.
2. Python programs
  - plot_meteors.py  --> Reads the positiions of of bolides from USG sample and makes a map with such events
    - Inputs:
    - map.dbf, map.shp --> Contains information to plot Mexico map
    - meteors_database.tab, USG_meteors_database.tab --> Tables with meteors position
    - output:
    - meteors_map.pdf  --> map with events position
  - energy_factor.py --> Compares energies of events which appears in samples GLM and USG and computes corrective factor via linear regression. Release a graph with the results
  - TEC_detrender.py --> Detrends a TEC series using a Golay-Savitsky filter. See Pradipta et al. 2018
  - TEC_series.py    --> Takes the detrended signals by TEC_detrender.py and plots the resultant time series
  - wavelet_power.py --> Calculates the wavelet tranform of detrended TEC series and plot wavelet power spectrum
3. Useful data
  - meteors_database.tab --> contains table of GLM sample
  - USG_meteors_database.tab --> contains table of USG sample
  - meteors_map.pdf --> map with events positions
4. Paper
  Has three versions: for elsevier journal, mdpi and arXiv. Contains the same abstract, sections and so on, only the format required by the journal varies.
  
