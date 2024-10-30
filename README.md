# 2010-2021_timeseries_analysis
Code for plotting and estimating trends for GOM long term trends paper
Data can be downloaded from: 
https://doi.org/10.5061/dryad.9zw3r22n4 (pending release by Dryad)

The main script for plotting these timeseries is plot_timeseries_weeks_nature_dryad.

All other scripts are called by this main script.

Detection labels are provided in the Labels folder. Each file contains a vector, zID with 3 columns
- Detection event timestamp as a matlab datenumber
- Detection label, a number, which can be associated with a species ID using the associated row of mySpID.
- Detection label score. A Softmax score providing some information on the likelihood on the quality of the label. 
Values near 1 indicate that the classifier scored this event as more likely to be correct. 
Labels with scores near 0 are less likely to be correct.

Daily recording effort information is provided in .mat files in the Effort folder, combined in the code with additional deployment information stored in the code folde.

Output is saved to the TimeSeries folder, organized by site.

Species codes include:
Gg: Grampus griseus
Kspp: Kogia spp.
Md: Mesoplodon densirostris
Me: Mesoplodon europaeus
Pm: Physeter macrocephalus
UD: Unidentified delphinid
UD_LF: Low frequency delphind
Zc: Ziphius cavirostris
SS: Snapping shrimp
Sonar: Echosounders
Ship: Ship noise


