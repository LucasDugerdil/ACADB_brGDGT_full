# BRT calibrations for Arid Central Asian for paleo brGDGT 

## Project description
This GitHub project is associated to the publication of *Boosted Regression Trees machine-learning method drastically improves the brGDGT-based climate reconstruction in drylands.* in *XXXX* accessible [here](https://www.researchgate.net/profile/Lucas-Dugerdil?ev=hdr_xprf)

This R script is the full script develloped in the publication.
It is usefull to verify the replicability of this study.
The script could be modified for calibrations in other study areas.
This script is available on zenodo and citable at 
[![DOI](https://zenodo.org/badge/952259754.svg)](https://doi.org/10.5281/zenodo.16658140)

## How to install/run the ACADB brGDGT calibrations?
1. Install [R](https://larmarange.github.io/analyse-R/installation-de-R-et-RStudio.html)
2. It is easier to use [Rstudio](https://posit.co/downloads/)
3. Download this GitHub repository from ZIP file (by clicing on the green button `<> Code` beyond. 
## To run the full code
	- Open the `ACADB_brGDGT_full.Rproj` file in Rstudio
	- (Optional) to apply the `randomTF()` test on your different models, turn to `TRUE` the test change into `test.randomTF = T`, then the script will launch the function `Plot.randomTF()`
	- (Optional) many option and settings for each function (e.g. change the `BRT` model for the `RF`, export figure in `plotly`, etc.) can be discovered when look at the `./Import/Script/BRT_script.R` file
