# GRAB5HT Sensor Imaging
This Repository contains the codes used to process and analyze in vivo fiber photometry GRAB5HT imaging dataset. Run in MATLAB.  

batch_unmixing.m	%% (function) Used to unmix raw data acquired from spectrometer. 

Indv_XXX.m	% (function) Used to pre-process unmixed data, including stitching data of baseline and treatment, applying polyfit to remove noise and artifacts, exporting plots for visualization.  	

GRAB5HT_XXX_var.m %% (function) Used to detect featured variables for each dataset.

GRAB5HT_XXX_bp.m %% (function) Used for batch processing data under certain folder.

GRAB5HT_excecute.m %% Master script to run analysis and export data. 

