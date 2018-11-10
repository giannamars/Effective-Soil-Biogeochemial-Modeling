# Effective-Soil-Biogeochemial-Modeling
Supplementary Material for "Equifinality, sloppiness, and emergent structures of soil biogeochemical models" by Marschmann et al. (2018)

* to reproduce *Fig.3* run PerfectDataSpec.py in Subfolders of Universally_Sloppy. Then call plot_spectra.py

* for *Fig.2* run the following scripts in Monod_Model: monod_model.py, generate_ensemble_predictions.py, generate_geodesics.py. Then plot with plot_timeseries.py, plot_landscape.py, plot_manifold_PCA.py

* for *Fig.4* run PlotSpectrum.py in the respective subfolders Functional_Gene_Data, Bulk_Data and Input_Output_Data. In order to check every iteration in the model reduction, check modelN_fit.py and modelN_reduce.py in the corresponding folder ModelN

* for Fig. 5 run Performance/plot_calibrations.py. For model predictions, dm me for the data files or uncomment lines in Functional_Gene_Data/Model0/model0_fit.py Functional_Gene_Data/Model32/model32_fit.py to run MCMC optimization ($/approx$ 2 days)
