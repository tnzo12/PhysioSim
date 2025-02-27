# PhysioSim
This physiologically based pharmacokinetic (PBPK) model is designed as an educational tool, hence named PBPKedu. The aim is to provide a simple, open source, freely downloadable PBPK model for performing basic single subject bottom-up simulations in R. The model has not been tested against a large set of compounds. All parameters were taken from the literature. The PBPK model was implemented by using the [rxode2](https://cran.r-project.org/web/packages/rxode2/index.html) package and non-compartmental analysis (NCA) was performed with the [PKNCA](https://cran.r-project.org/web/packages/PKNCA/index.html) package.

## Files
- `server.R`: Shiny server function, it contains the "logics" of the applications (e.g., what happens when clicking a button).
- `ui.R`: Shiny user interface function, it contains the visual aspect of the application.
- `main_pbpk_acat.R`: R script that can be used to run the model outside the Shiny app. This file is not used in `server.R` or `ui.R`.
- `functions/PBPK_model_rxode.R`: it contains the ordinary differential equations of the PBPK and compartmental absorption & transit (CAT) models, written in `RxODE` format.
- `functions/import_param.R`: it contains the functions used to import and elaborate the PBPK and CAT model parameters.
- `functions/functions_plot4.R`: it contains the functions used to generate the plots.
- `data/PBPK_parameters`: it contains `.xlsx` files of the physiological parameters for all the species supported by the PBPKedu app. All the references are included in the `.xlsx` files.
- `data/library_drugs`:it contains `.xlsx` files of the drug physicochemical parameters and of drugs pharmacokinetics. All the references are included in the `.xlsx` files and in `data/library_drugs/readme.txt`.
- `unified_tissue_comp`: contains tissue composition data for calculating partition coefficient
## How to install on a local machine

```
install.packages(c("shiny","shinyBS","shinyjs","bslib","readxl","writexl","rxode2","dplyr","ggplot2","RColorBrewer","gridExtra","PKNCA","shinybusy"))
```
