# BRT calibrations for Arid Central Asian for paleo brGDGT (application)

## Overview
This GitHub project is associated to the publication of "*Boosted Regression Trees machine-learning method drastically improves the brGDGT-based climate reconstruction in drylands.*" published in *Paleoceanography and Paleoclimatology* in 2025 (Dugerdil et al., 2025c).

**Author**: **Lucas Dugerdil**<sup>1,2</sup>

**Affiliations**:
1. Univ. Lyon, ENS de Lyon, Université Lyon 1, CNRS, UMR 5276 LGL-TPE, F-69364, Lyon, France1
2. Université de Montpellier, CNRS, IRD, EPHE, UMR 5554 ISEM, Montpellier, France

**ORCID**: [0000-0003-0266-564X](https://orcid.org/0000-0003-0266-564X)

**Funding**: ANR, Grant [ANR‐22‐CE27‐0018](https://anr.fr/Project-ANR-22-CE27-0018) (STEPABILITY), Sébastien Joannin

**Open Access**:

<table width="100%">
  <tr>
    <td width="33.33%" align="left" valign="middle"><strong>Research article</strong></td>
    <td width="33.33%" align="left" valign="middle"><strong>Published release</strong></td>
    <td width="33.33%" align="left" valign="middle"><strong>Data repository</strong></td>
  </tr>
  <tr>
    <td width="33.33%" align="left" valign="middle">

[![Static Badge](https://img.shields.io/badge/DOI-10.1029%2F2025PA005214-yellow)](https://doi.org/10.1029/2025PA005214)

</td>
    <td width="33.33%" align="left" valign="middle">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16679065.svg)](https://doi.org/10.5281/zenodo.16679065)

</td>
    <td width="33.33%" align="left" valign="middle">

[![Static Badge](https://img.shields.io/badge/DOI-10.1594%2FPANGAEA.983391-green)](https://doi.org/10.1594/PANGAEA.983391)

</td>
  </tr>
</table>

## Description 
This R script is the full script develloped in the publication.
It is usefull to verify the replicability of this study.
The script could be modified for calibrations in other study areas.
For application of the machine-learning reconstructions based on ACADB calibration set, please refer to the user-friendly GitHub repository [/ACADB_brGDGT_calibrations][https://github.com/LucasDugerdil/ACADB_brGDGT_calibrations]. 

## How to install/run the ACADB brGDGT calibrations?
1. Install [R](https://larmarange.github.io/analyse-R/installation-de-R-et-RStudio.html)
2. It is easier to use [Rstudio](https://posit.co/downloads/)
3. Download this GitHub repository from ZIP file (by clicing on the green button `<> Code` beyond. 
## To run the full code
	- Open the `ACADB_brGDGT_full.Rproj` file in Rstudio
	- (Optional) to apply the `randomTF()` test on your different models, turn to `TRUE` the test change into `test.randomTF = T`, then the script will launch the function `Plot.randomTF()`
	- (Optional) many option and settings for each function (e.g. change the `BRT` model for the `RF`, export figure in `plotly`, etc.) can be discovered when look at the `./Import/Script/BRT_script.R` file
