Arid_rock_weathering
========

[![Twitter Follow](https://img.shields.io/twitter/follow/espadrine.svg?style=social&label=Follow)](https://twitter.com/RoeyAngel)   [![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

Data, sequence analysis ploting for figures included in the paper: [New insights into the role of epilithic biological crust in arid rock weathering]().

Abstract
--------

Overview
--------
    ├── Data          # Primary data
    │   ├── desiccation_data.csv
    │   ├── Isotopes_data.csv
    │   ├── Rock_weathering_metadata_RA.csv
    │   ├── Rock_weathering_new2_otuTab.txt
    │   └── Rock_weathering_new2_silva.nrv119.taxonomy
    ├── LICENCE.txt          # Copyright information
    ├── packrat/          # R packages needed by the RMD file
    ├── README.md          # Overview of the repo
    ├── references.bib          # Bibtex formatted refereces cited in the RMD file
    ├── Results         # Result files generated while running the RMD file
    │   ├── Rock_weathering_Diversity.csv
    │   ├── Rock_weathering_Diversity.Rds
    │   ├── Rock_weathering_filt3_GMPR_Order_CI.csv
    │   ├── Rock_weathering_filt3_GMPR_Order_Pvals.csv
    │   ├── Rock_weathering_filt3_GMPR_Phylum_CI.csv
    │   ├── Rock_weathering_filt3_GMPR_Phylum_Pvals.csv
    │   ├── Rock_weathering_Richness.csv
    │   └── Rock_weathering_Richness.Rds
    ├── Rock_wathering_analysis_files.Rproj          # R project file
    ├── Rock_weathering_analysis4paper.html          # HTML output of the RMD file
    ├── Rock_weathering_analysis4paper.md          # MD output of the RMD file
    ├── Rock_weathering_analysis4paper.RMD            # Executable R markdown script
    ├── Rock_weathering_cache/          # Cache files generated while running the RMD file
    ├── Rock_weathering_figures/          # Figure files generated while running the RMD file
    └── Rock_weathering_process_sequences.sh          # Shell excecutable for processing the sequence data

Running the analysis
--------
The RMD file is best executed using [knitr](https://yihui.name/knitr/) on [RStudio](https://www.rstudio.com/). The packrat folder holds the specific package versions used for the analysis.
