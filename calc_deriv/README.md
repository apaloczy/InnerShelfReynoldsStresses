## Codes for calculating derived variables from raw moored ADCP data

This directory contains codes that were used to process the raw mooring data. Most of them produce one or more of the \*.npz or netCDF files located in the `../data_reproduce_figs` directory.

Brief description of each \*.m or \*.py file:

* `beamproc.m`: Loads raw ADCP data in along-beam coordinates and applies some QC steps.
* `enuproc.m`: Transforms the data saved by the `beamproc.m` script to velocities in Earth coordinates (east-north-up).
* `calculate_stresses.py`: Calculates Reynolds stresses from quality-controlled along-beam ADCP velocities, using the Dewey & Stringer (2020, unpublished manuscript) technique (see Equations 1 and 2 in the manuscript).
* `calc_noise_allADCPs.py`: Calculate noise floors for all five ADCPs (OC25SA, OC25SB, OC25M, OC40S and OC40N) from the spectra of along-beam velocities.

Please contact [André Palóczy](mailto:apaloczy@ucsd.edu) with any questions on the code.
