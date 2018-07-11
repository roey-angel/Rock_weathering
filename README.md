Arid_rock_weathering
==========

[![Twitter Follow](https://img.shields.io/twitter/follow/espadrine.svg?style=social&label=Follow)](https://twitter.com/RoeyAngel)

[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)  

Data, sequence analysis ploting for figures included in the paper: [New insights into the role of epilithic biological crust in arid rock weathering]().


**Abstract**

Overview
--------
    |- README          # Overview the repo
    |
    |- LICENSE         # Copyright information
    |
    |- .Rmd 	   # executable Rmarkdown for this study
    |
    |- .md 	   # 
    |
    |
    |- .csl 	   # Journal-specific formatting csl file
    |
    |- references.bib 	   # Bibtex formatted refereces cited in manuscript
    |
    |- data            # Raw and primary data
    |  |- README.md         # More specific overview of subdirectory
    |  |- metadata.tsv         # Metadata for mouse experimental groups
    |  |- kegg/  # Reference files used in analysis
    |  |- mapping/         # Normalized read counts mapping to *C. difficile* genes
    |  |- metabolic models/         # Output from bigSMALL analysis
    |  +- wetlab_assays/     # Raw data from wet lab experiments
    |
    |- code/           # Data analysis scripts
    |  |- README.md         # More specific overview of subdirectory
    |  |- R/           # R scipts for figures and analysis
    |      |- figures/      # Code to generatre main body figures
    |      |- supplement/      # Code to generatre supplementary figures
    |      |- support/      # Code to perform additional tasks needed for analysis
    |  |- python/      # python scripts
    |  +- pbs/         # pbs scripts
    |
    |- results         # All output from workflows and analyses
    |  |- README.md         # More specific overview of subdirectory
    |  |- tables/      # Text version of tables (main body)
    |  |- figures/     # Manuscript figures (main body)
    |  +- supplement/  # Supplementary data presentation
    |      |- figures/     # Supplementary figures
    |      +- tables/      # Excel versions of supplementary tables


#### Running the analysis
