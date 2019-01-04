---
title: "Role of BRC in arid rock weathering"
subtitle: "Data analysis and plotting for publication"
author: "Roey Angel (<roey.angel@bc.cas.cz>)"
date: "2019-01-04"
bibliography: references.bib
link-citations: yes
always_allow_html: yes
output:
  rmarkdown::html_document:
    toc: true
    toc_float: true
    keep_md: true
    keep_tex: true
    number_sections: false
    highlight: "pygments"
    theme: "flatly"
    dev: "png"
    df_print: "kable"
    fig_caption: true
    code_folding: "show"
---






## New insights into the role of epilithic biological crust in arid rock weathering
This script reproduces all sequence analysis steps and plots included in the paper plus some additional exploratory analyses.
The analysis is heavily based on the phyloseq package [@mcmurdie_phyloseq_2013], but also on many other R packages.


```r
set.seed(123456789)
bootstraps <- 1000
min_lib_size <- 1000
```

**Load data**

```r
read.csv("Data/Rock_weathering_new2_otuTab.txt", header = TRUE, row.names = 1, sep = "\t") %>%
  t() %>% 
  as.data.frame() ->
  Rock_weathering_OTUmat

sort_order <- as.numeric(gsub("OTU([0-9]+)", "\\1", colnames(Rock_weathering_OTUmat)))
Rock_weathering_OTUmat <- Rock_weathering_OTUmat[, order(sort_order)]
row.names(Rock_weathering_OTUmat) <- gsub("(.*)Nimrod[0-9]+|Osnat[0-9]+", "\\1", row.names(Rock_weathering_OTUmat))

Metadata <- read.csv("Data/Rock_weathering_metadata_RA.csv", row.names = 1, header = TRUE)
# Order abundance_mat samples according to the metadata
sample_order <- match(row.names(Rock_weathering_OTUmat), row.names(Metadata))
Rock_weathering_OTUmat %<>% arrange(., sample_order)
Metadata$sample_names <- row.names(Metadata)
Metadata$Uni.Source <- fct_collapse(Metadata$Source, Rock = c("Dolomite", "Limestone"))
Metadata$Climate.Source <-
  factor(
    paste(
      Metadata$Climate,
      Metadata$Source
    ),
    levels = c(
      "Arid Limestone",
      "Arid Dust",
      "Arid Loess soil",
      "Hyperarid Dolomite",
      "Hyperarid Dust",
      "Hyperarid Loess soil"
    ),
    labels = c(
      "Arid limestone",
      "Arid dust",
      "Arid loess soil",
      "Hyperarid dolomite",
      "Hyperarid dust",
      "Hyperarid loess soil"
    )
  )

Metadata$Climate.UniSource <-
  factor(
    paste(
      Metadata$Climate,
      Metadata$Uni.Source
    ),
    levels = c(
      "Arid Rock",
      "Arid Dust",
      "Arid Loess soil",
      "Hyperarid Rock",
      "Hyperarid Dust",
      "Hyperarid Loess soil"
    ),
    labels = c(
      "Arid rock",
      "Arid dust",
      "Arid loess soil",
      "Hyperarid rock",
      "Hyperarid dust",
      "Hyperarid loess soil"
    )
  )
# calculate sample size
Metadata$Lib.size = rowSums(Rock_weathering_OTUmat)
row.names(Rock_weathering_OTUmat) <- row.names(Metadata)

# Load taxonomy data
tax.file <- "Data/Rock_weathering_new2_silva.nrv119.taxonomy"
Taxonomy <- read.table(tax.file,  stringsAsFactors = FALSE) # read taxonomy file

# count how many ';' in each cell and add up to 6
for (i in 1:nrow(Taxonomy)) {
  semicolons <- length(gregexpr(";", Taxonomy$V2[i])[[1]])
  if (semicolons < 6) {
    x <- paste0(rep("Unclassified;", 6 - semicolons), collapse = "")
    Taxonomy$V2[i] <- paste0(Taxonomy$V2[i], x, sep = "")
  }
}

do.call( "rbind", strsplit( Taxonomy$V1, ";", fixed = TRUE)) %>% 
  gsub( "size=([0-9]+)", "\\1", .) %>%
  data.frame( ., do.call( "rbind", strsplit( Taxonomy$V2, ";", fixed = TRUE)), stringsAsFactors = F) %>% 
  apply(., 2, function(x) gsub( "\\(.*\\)", "", x)) %>% 
  replace(., . == "unclassified", "Unclassified") -> 
  Taxonomy

colnames( Taxonomy ) <- c( "OTU", "Frequency", "Domain", "Phylum", "Class", "Order", "Family", "Genus" )
# rownames(Taxonomy) <- colnames(Rock_weathering_OTUmat)
rownames(Taxonomy) <- Taxonomy[, 1]

# generate phyloseq object
Rock_dust <- phyloseq(otu_table(Rock_weathering_OTUmat, taxa_are_rows = FALSE),
                        tax_table(Taxonomy[, -c(1, 2)]),
                        sample_data(Metadata)
                        )

# Reorder factors for plotting
sample_data(Rock_dust)$Source %<>% fct_relevel("Limestone", "Dolomite", "Dust", "Loess soil")
```

Remove samples not for analysis

```r
samples2remove <- c(2, 3, 4, 5, 6, 7, 8, 10, 12) 
Rock_dust <- subset_samples(Rock_dust, !grepl(paste(c(sample_names(Rock_dust)[samples2remove]), collapse = "|"), sample_names(Rock_dust)))
Rock_dust <- filter_taxa(Rock_dust, function(x) sum(x) > 0, TRUE)

domains2remove <- c("", "Eukaryota", "Unclassified")
classes2remove <- c("Chloroplast")
families2remove <- c("Mitochondria")

Rock_weathering_filt <- subset_taxa(Rock_dust, !is.na(Phylum) &
                        !Domain %in% domains2remove &
                      !Class %in% classes2remove &
                      !Family %in% families2remove)
```

First let's explore the prevalence of different taxa in the database.

```r
prevdf <- apply(X = otu_table(Rock_weathering_filt),
                 MARGIN = ifelse(taxa_are_rows(Rock_weathering_filt), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(Rock_weathering_filt),
                      tax_table(Rock_weathering_filt))

prevdf %>%
  group_by(Phylum) %>%
  summarise(`Mean prevalence` = mean(Prevalence),
            `Sum prevalence` = sum(Prevalence)) ->
  Prevalence_phylum_summary

Prevalence_phylum_summary %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Phylum </th>
   <th style="text-align:right;"> Mean prevalence </th>
   <th style="text-align:right;"> Sum prevalence </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Acidobacteria </td>
   <td style="text-align:right;"> 16.9 </td>
   <td style="text-align:right;"> 1281 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Actinobacteria </td>
   <td style="text-align:right;"> 20.1 </td>
   <td style="text-align:right;"> 5660 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Aquificae </td>
   <td style="text-align:right;"> 14.3 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Armatimonadetes </td>
   <td style="text-align:right;"> 15.8 </td>
   <td style="text-align:right;"> 79 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidetes </td>
   <td style="text-align:right;"> 16.2 </td>
   <td style="text-align:right;"> 2337 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldiserica </td>
   <td style="text-align:right;"> 2.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Candidate_division_BRC1 </td>
   <td style="text-align:right;"> 14.0 </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Candidate_division_OD1 </td>
   <td style="text-align:right;"> 20.2 </td>
   <td style="text-align:right;"> 81 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Candidate_division_OP11 </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Candidate_division_TM7 </td>
   <td style="text-align:right;"> 18.3 </td>
   <td style="text-align:right;"> 440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlorobi </td>
   <td style="text-align:right;"> 15.5 </td>
   <td style="text-align:right;"> 62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chloroflexi </td>
   <td style="text-align:right;"> 14.9 </td>
   <td style="text-align:right;"> 2753 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cyanobacteria </td>
   <td style="text-align:right;"> 19.4 </td>
   <td style="text-align:right;"> 874 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Deinococcus-Thermus </td>
   <td style="text-align:right;"> 19.6 </td>
   <td style="text-align:right;"> 274 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elusimicrobia </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fibrobacteres </td>
   <td style="text-align:right;"> 16.2 </td>
   <td style="text-align:right;"> 65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Firmicutes </td>
   <td style="text-align:right;"> 16.2 </td>
   <td style="text-align:right;"> 1164 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fusobacteria </td>
   <td style="text-align:right;"> 14.0 </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatimonadetes </td>
   <td style="text-align:right;"> 15.4 </td>
   <td style="text-align:right;"> 848 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrospirae </td>
   <td style="text-align:right;"> 22.0 </td>
   <td style="text-align:right;"> 44 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NPL-UPA2 </td>
   <td style="text-align:right;"> 6.0 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Planctomycetes </td>
   <td style="text-align:right;"> 12.7 </td>
   <td style="text-align:right;"> 420 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proteobacteria </td>
   <td style="text-align:right;"> 19.1 </td>
   <td style="text-align:right;"> 4979 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SBYG-2791 </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SM2F11 </td>
   <td style="text-align:right;"> 20.0 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Spirochaetae </td>
   <td style="text-align:right;"> 15.0 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synergistetes </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tenericutes </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thaumarchaeota </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Verrucomicrobia </td>
   <td style="text-align:right;"> 13.2 </td>
   <td style="text-align:right;"> 383 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WCHB1-60 </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 57 </td>
  </tr>
</tbody>
</table>

```r
prevdf %>%
  group_by(Order) %>%
  summarise(`Mean prevalence` = mean(Prevalence),
            `Sum prevalence` = sum(Prevalence)) ->
  Prevalence_Order_summary

Prevalence_Order_summary %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Order </th>
   <th style="text-align:right;"> Mean prevalence </th>
   <th style="text-align:right;"> Sum prevalence </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 11B-2 </td>
   <td style="text-align:right;"> 9.2 </td>
   <td style="text-align:right;"> 55 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidimicrobiales </td>
   <td style="text-align:right;"> 18.2 </td>
   <td style="text-align:right;"> 709 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidithiobacillales </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Actinomycetales </td>
   <td style="text-align:right;"> 12.0 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Aeromonadales </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AKIW781 </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 592 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AKYG1722 </td>
   <td style="text-align:right;"> 17.4 </td>
   <td style="text-align:right;"> 157 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alteromonadales </td>
   <td style="text-align:right;"> 14.5 </td>
   <td style="text-align:right;"> 29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Anaerolineales </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 63 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Aquificales </td>
   <td style="text-align:right;"> 14.3 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ardenticatenales </td>
   <td style="text-align:right;"> 9.8 </td>
   <td style="text-align:right;"> 98 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AT425-EubC11_terrestrial_group </td>
   <td style="text-align:right;"> 18.2 </td>
   <td style="text-align:right;"> 328 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacillales </td>
   <td style="text-align:right;"> 18.3 </td>
   <td style="text-align:right;"> 475 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidales </td>
   <td style="text-align:right;"> 10.7 </td>
   <td style="text-align:right;"> 192 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BD2-11_terrestrial_group </td>
   <td style="text-align:right;"> 26.0 </td>
   <td style="text-align:right;"> 26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BD72BR169 </td>
   <td style="text-align:right;"> 14.0 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bdellovibrionales </td>
   <td style="text-align:right;"> 14.5 </td>
   <td style="text-align:right;"> 58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BG.g7 </td>
   <td style="text-align:right;"> 26.0 </td>
   <td style="text-align:right;"> 26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bifidobacteriales </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Brocadiales </td>
   <td style="text-align:right;"> 15.5 </td>
   <td style="text-align:right;"> 31 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Burkholderiales </td>
   <td style="text-align:right;"> 19.3 </td>
   <td style="text-align:right;"> 559 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C0119 </td>
   <td style="text-align:right;"> 12.0 </td>
   <td style="text-align:right;"> 120 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caenarcaniphilales </td>
   <td style="text-align:right;"> 23.0 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldilineales </td>
   <td style="text-align:right;"> 17.0 </td>
   <td style="text-align:right;"> 102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldisericales </td>
   <td style="text-align:right;"> 2.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Campylobacterales </td>
   <td style="text-align:right;"> 14.6 </td>
   <td style="text-align:right;"> 73 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caulobacterales </td>
   <td style="text-align:right;"> 23.4 </td>
   <td style="text-align:right;"> 257 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlorobiales </td>
   <td style="text-align:right;"> 15.5 </td>
   <td style="text-align:right;"> 62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chloroflexales </td>
   <td style="text-align:right;"> 7.4 </td>
   <td style="text-align:right;"> 52 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chromatiales </td>
   <td style="text-align:right;"> 16.8 </td>
   <td style="text-align:right;"> 84 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chthoniobacterales </td>
   <td style="text-align:right;"> 14.1 </td>
   <td style="text-align:right;"> 226 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Clostridiales </td>
   <td style="text-align:right;"> 15.4 </td>
   <td style="text-align:right;"> 401 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Corynebacteriales </td>
   <td style="text-align:right;"> 16.9 </td>
   <td style="text-align:right;"> 203 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:right;"> 18.6 </td>
   <td style="text-align:right;"> 1113 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Dehalococcoidales </td>
   <td style="text-align:right;"> 8.0 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Deinococcales </td>
   <td style="text-align:right;"> 19.6 </td>
   <td style="text-align:right;"> 255 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfobacterales </td>
   <td style="text-align:right;"> 8.4 </td>
   <td style="text-align:right;"> 42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfovibrionales </td>
   <td style="text-align:right;"> 8.0 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfurellales </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfuromonadales </td>
   <td style="text-align:right;"> 22.0 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elev-16S-976 </td>
   <td style="text-align:right;"> 23.5 </td>
   <td style="text-align:right;"> 47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EMP-G18 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Enterobacteriales </td>
   <td style="text-align:right;"> 27.0 </td>
   <td style="text-align:right;"> 135 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Erysipelotrichales </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 33 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Euzebyales </td>
   <td style="text-align:right;"> 19.4 </td>
   <td style="text-align:right;"> 252 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fibrobacterales </td>
   <td style="text-align:right;"> 16.2 </td>
   <td style="text-align:right;"> 65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Flavobacteriales </td>
   <td style="text-align:right;"> 16.4 </td>
   <td style="text-align:right;"> 164 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Frankiales </td>
   <td style="text-align:right;"> 28.2 </td>
   <td style="text-align:right;"> 367 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fusobacteriales </td>
   <td style="text-align:right;"> 14.0 </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gaiellales </td>
   <td style="text-align:right;"> 18.4 </td>
   <td style="text-align:right;"> 386 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gammaproteobacteria_Incertae_Sedis </td>
   <td style="text-align:right;"> 8.0 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatimonadales </td>
   <td style="text-align:right;"> 14.4 </td>
   <td style="text-align:right;"> 446 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GR-WP33-30 </td>
   <td style="text-align:right;"> 23.0 </td>
   <td style="text-align:right;"> 46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HOC36 </td>
   <td style="text-align:right;"> 8.0 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JG30-KF-CM45 </td>
   <td style="text-align:right;"> 20.8 </td>
   <td style="text-align:right;"> 604 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Kineosporiales </td>
   <td style="text-align:right;"> 24.2 </td>
   <td style="text-align:right;"> 145 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lactobacillales </td>
   <td style="text-align:right;"> 14.4 </td>
   <td style="text-align:right;"> 158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Legionellales </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lineage_IIb </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lineage_IV </td>
   <td style="text-align:right;"> 13.0 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LNR_A2-18 </td>
   <td style="text-align:right;"> 15.0 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Methylophilales </td>
   <td style="text-align:right;"> 11.5 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micrococcales </td>
   <td style="text-align:right;"> 24.3 </td>
   <td style="text-align:right;"> 365 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micromonosporales </td>
   <td style="text-align:right;"> 16.2 </td>
   <td style="text-align:right;"> 130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Myxococcales </td>
   <td style="text-align:right;"> 15.4 </td>
   <td style="text-align:right;"> 200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Neisseriales </td>
   <td style="text-align:right;"> 14.5 </td>
   <td style="text-align:right;"> 58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitriliruptorales </td>
   <td style="text-align:right;"> 21.3 </td>
   <td style="text-align:right;"> 64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrosomonadales </td>
   <td style="text-align:right;"> 28.0 </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrospirales </td>
   <td style="text-align:right;"> 22.0 </td>
   <td style="text-align:right;"> 44 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NKB5 </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Obscuribacterales </td>
   <td style="text-align:right;"> 12.5 </td>
   <td style="text-align:right;"> 50 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oceanospirillales </td>
   <td style="text-align:right;"> 15.0 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Opitutales </td>
   <td style="text-align:right;"> 10.0 </td>
   <td style="text-align:right;"> 30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Order_II </td>
   <td style="text-align:right;"> 16.0 </td>
   <td style="text-align:right;"> 32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Order_III </td>
   <td style="text-align:right;"> 16.0 </td>
   <td style="text-align:right;"> 48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Order_IV </td>
   <td style="text-align:right;"> 26.0 </td>
   <td style="text-align:right;"> 26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pasteurellales </td>
   <td style="text-align:right;"> 14.0 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Phycisphaerales </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Planctomycetales </td>
   <td style="text-align:right;"> 12.3 </td>
   <td style="text-align:right;"> 344 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Propionibacteriales </td>
   <td style="text-align:right;"> 19.6 </td>
   <td style="text-align:right;"> 352 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pseudomonadales </td>
   <td style="text-align:right;"> 19.6 </td>
   <td style="text-align:right;"> 372 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pseudonocardiales </td>
   <td style="text-align:right;"> 20.7 </td>
   <td style="text-align:right;"> 352 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PYR10d3 </td>
   <td style="text-align:right;"> 22.0 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhizobiales </td>
   <td style="text-align:right;"> 20.8 </td>
   <td style="text-align:right;"> 955 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhodobacterales </td>
   <td style="text-align:right;"> 21.9 </td>
   <td style="text-align:right;"> 219 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhodocyclales </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhodospirillales </td>
   <td style="text-align:right;"> 18.0 </td>
   <td style="text-align:right;"> 450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rickettsiales </td>
   <td style="text-align:right;"> 19.2 </td>
   <td style="text-align:right;"> 115 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rubrobacterales </td>
   <td style="text-align:right;"> 28.5 </td>
   <td style="text-align:right;"> 484 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S0134_terrestrial_group </td>
   <td style="text-align:right;"> 9.6 </td>
   <td style="text-align:right;"> 48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SC-I-84 </td>
   <td style="text-align:right;"> 15.0 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Selenomonadales </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Solirubrobacterales </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 989 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sphaerobacterales </td>
   <td style="text-align:right;"> 14.3 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sphingobacteriales </td>
   <td style="text-align:right;"> 15.3 </td>
   <td style="text-align:right;"> 751 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sphingomonadales </td>
   <td style="text-align:right;"> 24.0 </td>
   <td style="text-align:right;"> 552 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Streptomycetales </td>
   <td style="text-align:right;"> 23.0 </td>
   <td style="text-align:right;"> 69 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Streptosporangiales </td>
   <td style="text-align:right;"> 22.5 </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup_10 </td>
   <td style="text-align:right;"> 15.0 </td>
   <td style="text-align:right;"> 30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup_2 </td>
   <td style="text-align:right;"> 13.0 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup_3 </td>
   <td style="text-align:right;"> 16.1 </td>
   <td style="text-align:right;"> 129 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup_4 </td>
   <td style="text-align:right;"> 16.9 </td>
   <td style="text-align:right;"> 523 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup_6 </td>
   <td style="text-align:right;"> 19.2 </td>
   <td style="text-align:right;"> 441 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup_7 </td>
   <td style="text-align:right;"> 18.0 </td>
   <td style="text-align:right;"> 90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SubsectionI </td>
   <td style="text-align:right;"> 20.7 </td>
   <td style="text-align:right;"> 62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SubsectionII </td>
   <td style="text-align:right;"> 23.8 </td>
   <td style="text-align:right;"> 285 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SubsectionIII </td>
   <td style="text-align:right;"> 18.5 </td>
   <td style="text-align:right;"> 351 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SubsectionIV </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synergistales </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermales </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermoanaerobacterales </td>
   <td style="text-align:right;"> 13.2 </td>
   <td style="text-align:right;"> 53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermogemmatisporales </td>
   <td style="text-align:right;"> 14.5 </td>
   <td style="text-align:right;"> 29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermophilales </td>
   <td style="text-align:right;"> 22.8 </td>
   <td style="text-align:right;"> 91 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thiotrichales </td>
   <td style="text-align:right;"> 13.0 </td>
   <td style="text-align:right;"> 39 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TRA3-20 </td>
   <td style="text-align:right;"> 20.2 </td>
   <td style="text-align:right;"> 101 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Unclassified </td>
   <td style="text-align:right;"> 16.7 </td>
   <td style="text-align:right;"> 2355 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Unknown_Order </td>
   <td style="text-align:right;"> 14.6 </td>
   <td style="text-align:right;"> 73 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vampirovibrionales </td>
   <td style="text-align:right;"> 12.0 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Verrucomicrobiales </td>
   <td style="text-align:right;"> 8.0 </td>
   <td style="text-align:right;"> 24 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Xanthomonadales </td>
   <td style="text-align:right;"> 17.9 </td>
   <td style="text-align:right;"> 233 </td>
  </tr>
</tbody>
</table>

Based on that we'll remove all phyla with a prevalence of under 7

```r
Prevalence_phylum_summary %>% 
  filter(`Sum prevalence` < 7) %>% 
  select(Phylum) %>% 
  map(as.character) %>% 
  unlist() ->
  filterPhyla

Rock_weathering_filt2 <- subset_taxa(Rock_weathering_filt, !Phylum %in% filterPhyla)
sample_data(Rock_weathering_filt2)$Lib.size <- rowSums(otu_table(Rock_weathering_filt2))
print(Rock_weathering_filt)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1259 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 1259 taxa by 6 taxonomic ranks ]
```

```r
print(Rock_weathering_filt2)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1256 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 1256 taxa by 6 taxonomic ranks ]
```

Plot general prevalence features of the phyla

```r
# Subset to the remaining phyla
prevdf_phylum_filt <- subset(prevdf, Phylum %in% get_taxa_unique(Rock_weathering_filt2, "Phylum"))
ggplot(prevdf_phylum_filt,
       aes(TotalAbundance, Prevalence / nsamples(Rock_weathering_filt2), color = Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05,
             alpha = 0.5,
             linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Phylum) + theme(legend.position = "none")
```

![](Rock_weathering_figures/prevalence phylum-1.svg)<!-- -->

Plot general prevalence features of the top 20 orders

```r
# Subset to the remaining phyla
prevdf_order_filt <- subset(prevdf, Order %in% get_taxa_unique(Rock_weathering_filt2, "Order"))

# grab the top 30 most abundant orders
prevdf_order_filt %>% 
  group_by(Order) %>%
  summarise(Combined.abundance = sum(TotalAbundance)) %>% 
  arrange(desc(Combined.abundance)) %>% 
  .[1:30, "Order"]  ->
  Orders2plot

prevdf_order_filt2 <- subset(prevdf, Order %in% Orders2plot$Order)

ggplot(prevdf_order_filt2,
       aes(TotalAbundance, Prevalence / nsamples(Rock_weathering_filt2), color = Order)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05,
             alpha = 0.5,
             linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Order) + theme(legend.position = "none")
```

![](Rock_weathering_figures/prevalence order-1.svg)<!-- -->

#### Unsupervised filtering by prevalence
We'll remove all sequences which appear in less than 10% of the samples

```r
# Define prevalence threshold as 10% of total samples
prevalenceThreshold <- 0.1 * nsamples(Rock_weathering_filt)
prevalenceThreshold
```

```
## [1] 3.4
```

```r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <-
  row.names(prevdf_phylum_filt)[(prevdf_phylum_filt$Prevalence >= prevalenceThreshold)]
Rock_weathering_filt3 <- prune_taxa(keepTaxa, Rock_weathering_filt2)
sample_data(Rock_weathering_filt3)$Lib.size <- rowSums(otu_table(Rock_weathering_filt3))
print(Rock_weathering_filt2)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1256 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 1256 taxa by 6 taxonomic ranks ]
```

```r
print(Rock_weathering_filt3)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1249 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 1249 taxa by 6 taxonomic ranks ]
```
This removed 7 or 0.557% of the sequences.

### Exploring Rock_dust dataset features
First let's look at the count data distribution

```r
PlotLibDist(Rock_weathering_filt3)
```

![](Rock_weathering_figures/plot abundance-1.svg)<!-- -->

```r
sample_data(Rock_weathering_filt3) %>% 
  remove_rownames %>% 
  select(sample_title, Lib.size) %>% 
  as(., "data.frame") %>% 
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> sample_title </th>
   <th style="text-align:right;"> Lib.size </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Arid_Settled Dust_1 </td>
   <td style="text-align:right;"> 1562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Loess soil_8 </td>
   <td style="text-align:right;"> 6536 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Loess soil_10 </td>
   <td style="text-align:right;"> 3921 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Loess soil_12 </td>
   <td style="text-align:right;"> 4935 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Settled Dust_2 </td>
   <td style="text-align:right;"> 4421 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Settled Dust_1 </td>
   <td style="text-align:right;"> 1001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Settled Dust_2 </td>
   <td style="text-align:right;"> 16095 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_1 </td>
   <td style="text-align:right;"> 9765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_2 </td>
   <td style="text-align:right;"> 9130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_3 </td>
   <td style="text-align:right;"> 11218 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_4 </td>
   <td style="text-align:right;"> 13838 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_5 </td>
   <td style="text-align:right;"> 11177 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_6 </td>
   <td style="text-align:right;"> 10781 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_7 </td>
   <td style="text-align:right;"> 15417 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_8 </td>
   <td style="text-align:right;"> 9721 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_9 </td>
   <td style="text-align:right;"> 20927 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_10 </td>
   <td style="text-align:right;"> 16812 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_11 </td>
   <td style="text-align:right;"> 14325 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Limestone_12 </td>
   <td style="text-align:right;"> 5112 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_1 </td>
   <td style="text-align:right;"> 62166 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_2 </td>
   <td style="text-align:right;"> 73930 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_3 </td>
   <td style="text-align:right;"> 123438 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_4 </td>
   <td style="text-align:right;"> 74161 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_5 </td>
   <td style="text-align:right;"> 98998 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_6 </td>
   <td style="text-align:right;"> 97834 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_7 </td>
   <td style="text-align:right;"> 160207 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_8 </td>
   <td style="text-align:right;"> 78535 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_9 </td>
   <td style="text-align:right;"> 47155 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_10 </td>
   <td style="text-align:right;"> 52276 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_11 </td>
   <td style="text-align:right;"> 63267 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyperarid_Dolomite_12 </td>
   <td style="text-align:right;"> 53859 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Loess soil_1 </td>
   <td style="text-align:right;"> 61130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Loess soil_2 </td>
   <td style="text-align:right;"> 62204 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arid_Loess soil_3 </td>
   <td style="text-align:right;"> 55724 </td>
  </tr>
</tbody>
</table>
The figure and table indicate only a small deviation in the number of reads per samples.


```r
(mod1 <- adonis(
  otu_table(Rock_weathering_filt3) ~ Lib.size,
  data = as(sample_data(Rock_weathering_filt3), "data.frame"), 
  method = "bray",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3) ~ Lib.size,      data = as(sample_data(Rock_weathering_filt3), "data.frame"),      permutations = 9999, method = "bray") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Lib.size   1    2.5245 2.52447  8.7483 0.21469  1e-04 ***
## Residuals 32    9.2341 0.28857         0.78531           
## Total     33   11.7586                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
PlotReadHist(as(otu_table(Rock_weathering_filt3), "matrix"))
```

![](Rock_weathering_figures/mod abundance-1.svg)<!-- -->

```r
notAllZero <- (rowSums(t(otu_table(Rock_weathering_filt3))) > 0)
meanSdPlot(as.matrix(log2(t(otu_table(Rock_weathering_filt3))[notAllZero, ] + 1)))
```

![](Rock_weathering_figures/mod abundance-2.svg)<!-- -->

### Account for variation in library read-depth
We'll use the GMPR method [@chen_gmpr:_2017]

```r
Rock_weathering_filt3_GMPR <- Rock_weathering_filt3
Rock_weathering_filt3 %>%
  otu_table(.) %>%
  t() %>%
  as(., "matrix") %>%
  GMPR() ->
  GMPR_factors
```

```
## Begin GMPR size factor calculation ...
## Completed!
## Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers!
```

```r
Rock_weathering_filt3 %>%
  otu_table(.) %>%
  t() %*% diag(1 / GMPR_factors$gmpr) %>%
  t() %>%
  as.data.frame(., row.names = sample_names(Rock_weathering_filt3)) %>%
  otu_table(., taxa_are_rows = FALSE) ->
  otu_table(Rock_weathering_filt3_GMPR)
sample_data(Rock_weathering_filt3_GMPR)$Lib.size <- sample_sums(Rock_weathering_filt3_GMPR)

adonis(
  otu_table(Rock_weathering_filt3_GMPR) ~ Lib.size,
  data = as(sample_data(Rock_weathering_filt3_GMPR), "data.frame"),
  method = "bray",
  permutations = 9999
)
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR) ~ Lib.size,      data = as(sample_data(Rock_weathering_filt3_GMPR), "data.frame"),      permutations = 9999, method = "bray") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Lib.size   1    2.5834 2.58337  9.4557 0.22809  1e-04 ***
## Residuals 32    8.7426 0.27321         0.77191           
## Total     33   11.3260                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
PlotLibDist(Rock_weathering_filt3_GMPR)
```

![](Rock_weathering_figures/GMPR-1.svg)<!-- -->

```r
PlotReadHist(as(otu_table(Rock_weathering_filt3_GMPR), "matrix"))
```

![](Rock_weathering_figures/GMPR diag plots-1.svg)<!-- -->

```r
notAllZero <- (rowSums(t(otu_table(Rock_weathering_filt3_GMPR))) > 0)
meanSdPlot(as.matrix(log2(t(otu_table(Rock_weathering_filt3_GMPR))[notAllZero, ] + 1)))
```

![](Rock_weathering_figures/GMPR diag plots-2.svg)<!-- -->

### Alpha diversity 
Calculate and plot alpha diversity mertrics.

```r
# non-parametric richness estimates
rarefaction.mat <- matrix(0, nrow = nsamples(Rock_weathering_filt3), ncol = bootstraps)
rownames(rarefaction.mat) <- sample_names(Rock_weathering_filt3)
rich.ests <- list(S.obs = rarefaction.mat, S.chao1 = rarefaction.mat, se.chao1 = rarefaction.mat,
                   S.ACE = rarefaction.mat, se.ACE = rarefaction.mat)

for (i in seq(bootstraps)) {
  sub.OTUmat <- rrarefy(otu_table(Rock_weathering_filt3), min(rowSums(otu_table(Rock_weathering_filt3))))
  for (j in seq(length(rich.ests))) {
    rich.ests[[j]][, i] <- t(estimateR(sub.OTUmat))[, j]
  }
}

Richness <- data.frame(row.names = row.names(rich.ests[[1]]))
for (i in c(1, seq(2, length(rich.ests), 2))) {
  S <- apply(rich.ests[[i]], 1, mean)
  if (i == 1) { 
    se <- apply(rich.ests[[i]], 1, function(x) (mean(x)/sqrt(length(x))))
    } else se <- apply(rich.ests[[i + 1]], 1, mean)
  Richness <- cbind(Richness, S, se)
}
colnames(Richness) <- c("S.obs", "S.obs.se", "S.chao1", "S.chao1.se", "S.ACE", "S.ACE.se")


saveRDS(Richness, file = "Results/Rock_weathering_Richness.Rds")
write.csv(Richness, file = "Results/Rock_weathering_Richness.csv")

ses <- grep("\\.se", colnames(Richness))
Richness[, ses] %>% 
  gather(key = "est.se") -> se.dat
Richness[, -unique(ses)] %>% 
  gather(key = "est") -> mean.dat

n <- length(unique(mean.dat$est))

# diversity indices
diversity.inds <- list(Shannon = rarefaction.mat, inv.simpson = rarefaction.mat, BP = rarefaction.mat)
for (i in seq(bootstraps)) {
  sub.OTUmat <- rrarefy(otu_table(Rock_weathering_filt3), min(rowSums(otu_table(Rock_weathering_filt3))))
  diversity.inds$Shannon[, i] <- diversityresult(sub.OTUmat, index = 'Shannon', method = 'each site', digits = 3)[, 1]
  diversity.inds$inv.simpson[, i] <- diversityresult(sub.OTUmat, index = 'inverseSimpson', method = 'each site', digits = 3)[, 1]
  diversity.inds$BP[, i] <- diversityresult(sub.OTUmat, index = 'Berger', method = 'each site', digits = 3)[, 1]
}

Diversity <- data.frame(row.names = row.names(diversity.inds[[1]]))
for (i in seq(length(diversity.inds))) {
  S <- apply(diversity.inds[[i]], 1, mean)
  se <- apply(diversity.inds[[i]], 1, function(x) (mean(x)/sqrt(length(x))))
  Diversity <- cbind(Diversity, S, se)
}
colnames(Diversity) <- c("Shannon", "Shannon.se", "Inv.simpson", "Inv.simpson.se", "BP", "BP.se")

ses <- grep("\\.se", colnames(Diversity))
Diversity[, ses] %>% gather(key = "est.se") -> se.dat
Diversity[, -unique(ses)] %>% gather(key = "est") -> mean.dat

saveRDS(Diversity, file = "Results/Rock_weathering_Diversity.Rds")
write.csv(Diversity, file = "Results/Rock_weathering_Diversity.csv")
```

Test the differences in alpha diversity.


#### Plot all alpha diversity metrics together

```r
Richness_Diversity_long[Richness_Diversity_long$Metric != "Chao1" &
                          Richness_Diversity_long$Metric != "Inv. Simpson" &
                          Richness_Diversity_long$Metric != "Berger Parker", ] %>% 
  droplevels() ->
  Richness_Diversity_long2plot

p_alpha <- ggplot(Richness_Diversity_long2plot, aes(
  x = Source,
  y = Estimate
)) +
  geom_violin(aes(colour = Climate, fill = Climate), alpha = 1/3) +
  geom_jitter(aes(colour = Climate, fill = Climate), shape = 16, size = 2, width = 0.2, alpha = 2/3) +
  scale_colour_manual(values = pom4, name = "") +
  scale_fill_manual(values = pom4, name = "") +
  theme_cowplot(font_size = 11, font_family = f_name) +
  # geom_errorbar(alpha = 1 / 2, width = 0.3) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.9,
    hjust = 0.9
  )) +
  facet_grid(Metric ~ Climate, shrink = FALSE, scale = "free") +
  background_grid(major = "y",
                  minor = "none") +
  theme(panel.spacing = unit(2, "lines"))

dat_text <- data.frame(
  label = as.character(fct_c(ph_Sobs$groups$groups, ph_ACE$groups$groups, ph_Shannon$groups$groups)),
  Metric   = rep(levels(Richness_Diversity_long2plot$Metric), each = 6),
  Climate = str_split(rownames(ph_Sobs$groups), ":", simplify = TRUE)[, 1], 
  x = c("Loess soil", "Loess soil", "Limestone", "Dust", "Dolomite", "Dust"),
  # x     = as.factor(levels(Richness_Diversity_long2plot$Climate.Source)),
  y = rep(c(460, 850, 6.5), each = 6)
  # y = rep(c(40, 140, 0.5), each = 6)
)


p_alpha <- p_alpha + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label),
  nudge_x = -0.2,
  nudge_y = -0.1
)

print(p_alpha)
```

![](Rock_weathering_figures/plot alpha-1.svg)<!-- -->

```r
Richness_Diversity_long2plot %>%
  group_by(Metric, Climate.Source) %>%   # the grouping variable
  summarise(mean_PL = mean(Estimate),  # calculates the mean of each group
            sd_PL = sd(Estimate), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(Estimate)/sqrt(n())) %>% # calculates the standard error of each group
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Metric </th>
   <th style="text-align:left;"> Climate.Source </th>
   <th style="text-align:right;"> mean_PL </th>
   <th style="text-align:right;"> sd_PL </th>
   <th style="text-align:right;"> n_PL </th>
   <th style="text-align:right;"> SE_PL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S obs. </td>
   <td style="text-align:left;"> Arid limestone </td>
   <td style="text-align:right;"> 181.614750 </td>
   <td style="text-align:right;"> 51.1677757 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 14.7708645 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S obs. </td>
   <td style="text-align:left;"> Arid dust </td>
   <td style="text-align:right;"> 169.249500 </td>
   <td style="text-align:right;"> 109.0153596 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 77.0855000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S obs. </td>
   <td style="text-align:left;"> Arid loess soil </td>
   <td style="text-align:right;"> 416.015667 </td>
   <td style="text-align:right;"> 8.2772347 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 4.7788637 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S obs. </td>
   <td style="text-align:left;"> Hyperarid dolomite </td>
   <td style="text-align:right;"> 128.760500 </td>
   <td style="text-align:right;"> 31.2564362 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 9.0229559 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S obs. </td>
   <td style="text-align:left;"> Hyperarid dust </td>
   <td style="text-align:right;"> 107.261000 </td>
   <td style="text-align:right;"> 88.7263447 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 62.7390000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S obs. </td>
   <td style="text-align:left;"> Hyperarid loess soil </td>
   <td style="text-align:right;"> 220.405000 </td>
   <td style="text-align:right;"> 54.3993291 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 31.4074673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACE </td>
   <td style="text-align:left;"> Arid limestone </td>
   <td style="text-align:right;"> 353.781788 </td>
   <td style="text-align:right;"> 68.8104346 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 19.8638615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACE </td>
   <td style="text-align:left;"> Arid dust </td>
   <td style="text-align:right;"> 334.651440 </td>
   <td style="text-align:right;"> 143.6013458 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 101.5414854 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACE </td>
   <td style="text-align:left;"> Arid loess soil </td>
   <td style="text-align:right;"> 746.497431 </td>
   <td style="text-align:right;"> 20.7936696 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 12.0052307 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACE </td>
   <td style="text-align:left;"> Hyperarid dolomite </td>
   <td style="text-align:right;"> 314.743324 </td>
   <td style="text-align:right;"> 100.4297817 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 28.9915808 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACE </td>
   <td style="text-align:left;"> Hyperarid dust </td>
   <td style="text-align:right;"> 311.050796 </td>
   <td style="text-align:right;"> 182.7650780 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 129.2344260 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACE </td>
   <td style="text-align:left;"> Hyperarid loess soil </td>
   <td style="text-align:right;"> 466.260376 </td>
   <td style="text-align:right;"> 42.5794120 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 24.5832350 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Shannon </td>
   <td style="text-align:left;"> Arid limestone </td>
   <td style="text-align:right;"> 3.782559 </td>
   <td style="text-align:right;"> 1.0939745 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3158033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Shannon </td>
   <td style="text-align:left;"> Arid dust </td>
   <td style="text-align:right;"> 2.997964 </td>
   <td style="text-align:right;"> 1.5634435 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1.1055215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Shannon </td>
   <td style="text-align:left;"> Arid loess soil </td>
   <td style="text-align:right;"> 5.594932 </td>
   <td style="text-align:right;"> 0.0319481 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.0184452 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Shannon </td>
   <td style="text-align:left;"> Hyperarid dolomite </td>
   <td style="text-align:right;"> 3.328207 </td>
   <td style="text-align:right;"> 0.2683934 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.0774785 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Shannon </td>
   <td style="text-align:left;"> Hyperarid dust </td>
   <td style="text-align:right;"> 1.483307 </td>
   <td style="text-align:right;"> 0.8891069 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.6286935 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Shannon </td>
   <td style="text-align:left;"> Hyperarid loess soil </td>
   <td style="text-align:right;"> 3.779445 </td>
   <td style="text-align:right;"> 0.8513243 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.4915123 </td>
  </tr>
</tbody>
</table>

### Beta diversity
Calculate and plot beta diversity mertrics. 

#### 1. Possible changes in biofilm arid vs hyper arid
Is there a difference between the two sites. However, since we know that that samples are of different nature we'll have to control for rock type, source and location:


```r
(mod1 <-  adonis(
  otu_table(Rock_weathering_filt3_GMPR) ~ Climate * Source * Location,
  data = as(sample_data(Rock_weathering_filt3_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR) ~ Climate *      Source * Location, data = as(sample_data(Rock_weathering_filt3_GMPR),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Climate           1    2.5299 2.52988 21.3820 0.21632 0.0001 ***
## Source            3    4.7060 1.56865 13.2579 0.40238 0.0001 ***
## Location          1    0.4870 0.48696  4.1157 0.04164 0.0033 ** 
## Climate:Source    1    0.4397 0.43972  3.7164 0.03760 0.0026 ** 
## Climate:Location  1    0.4605 0.46052  3.8922 0.03938 0.0056 ** 
## Source:Location   1    0.1142 0.11421  0.9653 0.00977 0.4812    
## Residuals        25    2.9580 0.11832         0.25292           
## Total            33   11.6952                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
Rock_weathering_filt3_GMPR_Arid <- subset_samples(Rock_weathering_filt3_GMPR, Climate == "Arid")
Rock_weathering_filt3_GMPR_Arid <- filter_taxa(Rock_weathering_filt3_GMPR_Arid, function(x) sum(x) > 0, TRUE)
(mod2 <-  adonis(
  otu_table(Rock_weathering_filt3_GMPR_Arid) ~ Source * Location,
  data = as(sample_data(Rock_weathering_filt3_GMPR_Arid), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR_Arid) ~      Source * Location, data = as(sample_data(Rock_weathering_filt3_GMPR_Arid),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Source     2    2.3058 1.15288  7.3438 0.47881 0.0001 ***
## Location   1    0.4691 0.46905  2.9878 0.09740 0.0132 *  
## Residuals 13    2.0408 0.15699         0.42379           
## Total     16    4.8156                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
Rock_weathering_filt3_GMPR_Hyperarid <- subset_samples(Rock_weathering_filt3_GMPR, Climate == "Hyperarid")
Rock_weathering_filt3_GMPR_Hyperarid <- filter_taxa(Rock_weathering_filt3_GMPR_Hyperarid, function(x) sum(x) > 0, TRUE)
(mod3 <-  adonis(
  otu_table(Rock_weathering_filt3_GMPR_Hyperarid) ~ Source * Location,
  data = as(sample_data(Rock_weathering_filt3_GMPR_Hyperarid), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR_Hyperarid) ~      Source * Location, data = as(sample_data(Rock_weathering_filt3_GMPR_Hyperarid),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##                 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Source           2    2.8454 1.42270 18.6150 0.65416 0.0001 ***
## Location         1    0.4729 0.47295  6.1882 0.10873 0.0102 *  
## Source:Location  1    0.1142 0.11421  1.4944 0.02626 0.2334    
## Residuals       12    0.9171 0.07643         0.21085           
## Total           16    4.3497                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
According to this model we see that indeed there's an effect of site on the community (p = 0.001), and that effect accounts for about 17% of the variance. Also, considering that Location is only borderline significant and explains very little of the data, we could probably take it out of the model to make a minimal adequate model.


```r
(mod4 <-  adonis(
  otu_table(Rock_weathering_filt3_GMPR) ~ Climate * Source,
  data = as(sample_data(Rock_weathering_filt3_GMPR), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR) ~ Climate *      Source, data = as(sample_data(Rock_weathering_filt3_GMPR),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Climate         1    2.5299 2.52988 17.6466 0.21632  1e-04 ***
## Source          3    4.7060 1.56865 10.9418 0.40238  1e-04 ***
## Climate:Source  1    0.4452 0.44521  3.1055 0.03807  8e-03 ** 
## Residuals      28    4.0142 0.14336         0.34323           
## Total          33   11.6952                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(mod5 <-  adonis(
  otu_table(Rock_weathering_filt3_GMPR_Arid) ~ Source,
  data = as(sample_data(Rock_weathering_filt3_GMPR_Arid), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR_Arid) ~      Source, data = as(sample_data(Rock_weathering_filt3_GMPR_Arid),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Source     2    2.3058 1.15288  6.4307 0.47881  2e-04 ***
## Residuals 14    2.5099 0.17928         0.52119           
## Total     16    4.8156                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(mod6 <-  adonis(
  otu_table(Rock_weathering_filt3_GMPR_Hyperarid) ~ Source,
  data = as(sample_data(Rock_weathering_filt3_GMPR_Hyperarid), "data.frame"),
  method = "horn",
  permutations = 9999
))
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR_Hyperarid) ~      Source, data = as(sample_data(Rock_weathering_filt3_GMPR_Hyperarid),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Source     2    2.8454 1.42270  13.241 0.65416  2e-04 ***
## Residuals 14    1.5043 0.10745         0.34584           
## Total     16    4.3497                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Final model

```r
print(mod4)
```

```
## 
## Call:
## adonis(formula = otu_table(Rock_weathering_filt3_GMPR) ~ Climate *      Source, data = as(sample_data(Rock_weathering_filt3_GMPR),      "data.frame"), permutations = 9999, method = "horn") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Climate         1    2.5299 2.52988 17.6466 0.21632  1e-04 ***
## Source          3    4.7060 1.56865 10.9418 0.40238  1e-04 ***
## Climate:Source  1    0.4452 0.44521  3.1055 0.03807  8e-03 ** 
## Residuals      28    4.0142 0.14336         0.34323           
## Total          33   11.6952                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
mod4_pairwise <- PairwiseAdonis(
  otu_table(Rock_weathering_filt3_GMPR),
  sample_data(Rock_weathering_filt3_GMPR)$Climate.Source,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod4_pairwise)
```

```
##                                         pairs total.DF    F.Model        R2   p.value
## 1           Arid dust vs Hyperarid loess soil        4   4.803200 0.6155424 0.1000000
## 2                 Arid dust vs Hyperarid dust        3   2.992098 0.5993669 0.3333333
## 3                 Arid dust vs Arid limestone       13   4.606432 0.2773884 0.0117000
## 4             Arid dust vs Hyperarid dolomite       13   6.594183 0.3546369 0.0119000
## 5                Arid dust vs Arid loess soil        4   5.363663 0.6413055 0.1000000
## 6      Hyperarid loess soil vs Hyperarid dust        4  10.866238 0.7836472 0.2000000
## 7      Hyperarid loess soil vs Arid limestone       14   5.593760 0.3008407 0.0070000
## 8  Hyperarid loess soil vs Hyperarid dolomite       14  16.022459 0.5520710 0.0026000
## 9     Hyperarid loess soil vs Arid loess soil        5  95.510242 0.9598031 0.1000000
## 10           Hyperarid dust vs Arid limestone       13   5.343419 0.3080949 0.0313000
## 11       Hyperarid dust vs Hyperarid dolomite       13  11.619149 0.4919377 0.0094000
## 12          Hyperarid dust vs Arid loess soil        4 205.847549 0.9856355 0.1000000
## 13       Arid limestone vs Hyperarid dolomite       23  21.900583 0.4988677 0.0001000
## 14          Arid limestone vs Arid loess soil       14   9.112582 0.4120994 0.0018000
## 15      Hyperarid dolomite vs Arid loess soil       14  16.841743 0.5643686 0.0032000
##    p.adjusted sig
## 1  0.11538462    
## 2  0.33333333    
## 3  0.02231250   .
## 4  0.02231250   .
## 5  0.11538462    
## 6  0.21428571    
## 7  0.02100000   .
## 8  0.01200000   .
## 9  0.11538462    
## 10 0.05216667    
## 11 0.02231250   .
## 12 0.11538462    
## 13 0.00150000   *
## 14 0.01200000   .
## 15 0.01200000   .
```

```r
sig_pairs <- as.character(mod4_pairwise$pairs[mod4_pairwise$p.adjusted < 0.05])
  
simper(otu_table(Rock_weathering_filt3_GMPR), sample_data(Rock_weathering_filt3_GMPR)$Climate.Source, parallel = 4)
```

```
## cumulative contributions of most influential species:
## 
## $`Arid dust_Hyperarid loess soil`
##      OTU6     OTU65    OTU838     OTU90    OTU596     OTU11    OTU187     OTU93    OTU746 
## 0.2223516 0.3677866 0.4056901 0.4304991 0.4528999 0.4693262 0.4837166 0.4963321 0.5089143 
##    OTU167    OTU711    OTU121    OTU144     OTU99    OTU194    OTU105    OTU356    OTU340 
## 0.5206772 0.5323162 0.5428037 0.5526841 0.5624404 0.5714645 0.5787230 0.5859611 0.5926902 
##    OTU115    OTU640     OTU55    OTU298    OTU715     OTU88    OTU197     OTU48    OTU221 
## 0.5992910 0.6058761 0.6123625 0.6187386 0.6250711 0.6313974 0.6376038 0.6433882 0.6488153 
##    OTU301   OTU1047    OTU322     OTU16    OTU586     OTU46    OTU172    OTU386     OTU67 
## 0.6541319 0.6593526 0.6645414 0.6694651 0.6739228 0.6782800 0.6825173 0.6866930 0.6908583 
##    OTU333    OTU854    OTU883 
## 0.6945572 0.6981826 0.7016897 
## 
## $`Arid dust_Hyperarid dust`
##      OTU6     OTU65    OTU838     OTU55 
## 0.4922822 0.6517818 0.6971512 0.7299602 
## 
## $`Arid dust_Arid limestone`
##      OTU6     OTU65     OTU20     OTU40    OTU838    OTU225    OTU422    OTU596     OTU73 
## 0.1123216 0.2223095 0.2918201 0.3566264 0.3846318 0.4019503 0.4188998 0.4354418 0.4510851 
##     OTU46     OTU11    OTU388    OTU119     OTU29     OTU18     OTU34    OTU746    OTU711 
## 0.4642743 0.4766768 0.4880677 0.4986328 0.5091612 0.5193871 0.5289518 0.5383129 0.5469937 
##     OTU16    OTU235    OTU854    OTU174    OTU299    OTU194     OTU68    OTU545     OTU37 
## 0.5556579 0.5642893 0.5727110 0.5809073 0.5890471 0.5960083 0.6028022 0.6094612 0.6160857 
##    OTU925     OTU41    OTU115    OTU140    OTU356     OTU60    OTU137    OTU197    OTU166 
## 0.6225697 0.6287094 0.6345215 0.6403150 0.6457386 0.6510360 0.6562074 0.6613174 0.6663109 
##    OTU605     OTU69     OTU55   OTU1089     OTU27    OTU218     OTU48     OTU45 
## 0.6712870 0.6762240 0.6810537 0.6854239 0.6895759 0.6936273 0.6976757 0.7016868 
## 
## $`Arid dust_Hyperarid dolomite`
##      OTU1      OTU2      OTU3      OTU7      OTU5     OTU65      OTU9     OTU18      OTU8 
## 0.1527495 0.2363471 0.3165709 0.3820169 0.4147143 0.4454892 0.4704565 0.4947935 0.5170144 
##     OTU17     OTU41     OTU16     OTU11    OTU936     OTU20     OTU28     OTU12     OTU10 
## 0.5378908 0.5586779 0.5762334 0.5933941 0.6090247 0.6239327 0.6383626 0.6526239 0.6663706 
##      OTU6     OTU15     OTU33     OTU14 
## 0.6767124 0.6869395 0.6964878 0.7058854 
## 
## $`Arid dust_Arid loess soil`
##      OTU65      OTU25      OTU88      OTU62      OTU11     OTU838      OTU68     OTU596 
## 0.07600412 0.09862665 0.12048117 0.14055785 0.15967703 0.17826462 0.19067874 0.20145631 
##     OTU177      OTU78     OTU133     OTU144      OTU46     OTU115      OTU67      OTU16 
## 0.21192993 0.22239646 0.23195007 0.24115456 0.25014254 0.25870758 0.26685638 0.27420872 
##     OTU156     OTU388     OTU116     OTU137     OTU197     OTU100      OTU91     OTU687 
## 0.28151396 0.28845987 0.29497148 0.30133776 0.30763580 0.31387308 0.31978477 0.32560209 
##     OTU746     OTU711      OTU73     OTU412     OTU194     OTU422     OTU160     OTU218 
## 0.33120586 0.33673375 0.34225822 0.34731488 0.35212805 0.35690940 0.36151216 0.36607842 
##     OTU152      OTU53     OTU131     OTU135     OTU356     OTU125     OTU163      OTU37 
## 0.37046694 0.37485116 0.37885223 0.38274887 0.38655530 0.39035864 0.39395340 0.39753331 
##     OTU556     OTU311     OTU130     OTU214     OTU227     OTU240     OTU341     OTU233 
## 0.40109574 0.40463803 0.40814028 0.41156987 0.41498306 0.41836142 0.42169721 0.42499571 
##     OTU854     OTU648     OTU101      OTU92      OTU55     OTU716     OTU220     OTU139 
## 0.42827522 0.43153536 0.43476556 0.43797394 0.44115558 0.44428602 0.44740514 0.45049071 
##     OTU140     OTU278      OTU20     OTU103      OTU60     OTU253     OTU171    OTU1235 
## 0.45356227 0.45649184 0.45941847 0.46234298 0.46512485 0.46790496 0.47063934 0.47334811 
##     OTU107     OTU365     OTU294     OTU301     OTU647     OTU206     OTU136     OTU231 
## 0.47601073 0.47866935 0.48128416 0.48386623 0.48644557 0.48901761 0.49157714 0.49411399 
##     OTU479    OTU1050    OTU1288      OTU82     OTU248     OTU360     OTU181     OTU434 
## 0.49659608 0.49905562 0.50151265 0.50394266 0.50636675 0.50871712 0.51105399 0.51337992 
##     OTU945     OTU638    OTU1102     OTU315    OTU1184     OTU457      OTU44     OTU491 
## 0.51570002 0.51801894 0.52031981 0.52260630 0.52487989 0.52712302 0.52935591 0.53158523 
##     OTU472     OTU938     OTU184     OTU471     OTU182     OTU538     OTU352     OTU169 
## 0.53380813 0.53602510 0.53823459 0.54042578 0.54253898 0.54464035 0.54673555 0.54882651 
##     OTU586     OTU191     OTU190     OTU751     OTU304     OTU251     OTU164     OTU667 
## 0.55091733 0.55300683 0.55508810 0.55712859 0.55915785 0.56117530 0.56318056 0.56518015 
##      OTU69     OTU619     OTU170     OTU154     OTU992    OTU1013     OTU707     OTU302 
## 0.56717714 0.56916832 0.57115385 0.57313482 0.57509789 0.57705219 0.57900145 0.58092857 
##     OTU725     OTU430    OTU1116     OTU132     OTU466     OTU364     OTU489     OTU166 
## 0.58282302 0.58471488 0.58660087 0.58848567 0.59036511 0.59223482 0.59409227 0.59593839 
##      OTU33    OTU1218     OTU421    OTU1319     OTU545     OTU309     OTU456     OTU562 
## 0.59778312 0.59962319 0.60146290 0.60326476 0.60505799 0.60683960 0.60859329 0.61034045 
##     OTU635     OTU468     OTU610     OTU347     OTU833     OTU582      OTU41     OTU333 
## 0.61206287 0.61376517 0.61546266 0.61714945 0.61882254 0.62049477 0.62216179 0.62381474 
##       OTU6      OTU27      OTU86     OTU451     OTU162     OTU561     OTU439     OTU579 
## 0.62545802 0.62705709 0.62865264 0.63024097 0.63181778 0.63338508 0.63491568 0.63644400 
##     OTU255     OTU653     OTU279     OTU368     OTU783      OTU93     OTU298     OTU323 
## 0.63797113 0.63949107 0.64100660 0.64251744 0.64400822 0.64549236 0.64697241 0.64843615 
##    OTU1024    OTU1005     OTU389     OTU799     OTU348     OTU543     OTU524     OTU258 
## 0.64988797 0.65133108 0.65277123 0.65418616 0.65559921 0.65699846 0.65838960 0.65977752 
##     OTU953     OTU929     OTU174    OTU1100     OTU999     OTU705     OTU199    OTU1049 
## 0.66116274 0.66251874 0.66386908 0.66520644 0.66653398 0.66783135 0.66912435 0.67041127 
##    OTU1242     OTU449     OTU303     OTU232     OTU765     OTU499      OTU10    OTU1011 
## 0.67169466 0.67295718 0.67420516 0.67544956 0.67668808 0.67792553 0.67913601 0.68034452 
##     OTU595     OTU697     OTU548     OTU353     OTU310      OTU28    OTU1170     OTU404 
## 0.68154983 0.68274660 0.68394090 0.68513515 0.68632833 0.68750914 0.68868372 0.68984795 
##     OTU141     OTU277     OTU469     OTU883     OTU557      OTU95      OTU29    OTU1108 
## 0.69101177 0.69216462 0.69330654 0.69444839 0.69558676 0.69671804 0.69784612 0.69896224 
##     OTU641 
## 0.70007721 
## 
## $`Hyperarid loess soil_Hyperarid dust`
##      OTU6     OTU55     OTU90     OTU11     OTU94    OTU187     OTU93    OTU167    OTU121 
## 0.4547139 0.4932104 0.5179318 0.5361943 0.5506355 0.5649686 0.5775294 0.5893159 0.5997621 
##     OTU99    OTU144      OTU9    OTU105    OTU340    OTU298    OTU115    OTU640    OTU715 
## 0.6095933 0.6193616 0.6270627 0.6343667 0.6410613 0.6476927 0.6542631 0.6608232 0.6672825 
##     OTU88    OTU197    OTU221    OTU322   OTU1047     OTU16    OTU484 
## 0.6735780 0.6797532 0.6851539 0.6903247 0.6951534 0.6999152 0.7044646 
## 
## $`Hyperarid loess soil_Arid limestone`
##      OTU6     OTU20     OTU40     OTU90    OTU225     OTU11     OTU73    OTU422    OTU119 
## 0.2070926 0.2740292 0.3382045 0.3562460 0.3732290 0.3885619 0.4034598 0.4178979 0.4286859 
##     OTU29    OTU187     OTU46    OTU388     OTU18    OTU167     OTU34     OTU93    OTU299 
## 0.4392449 0.4494930 0.4594962 0.4689023 0.4780078 0.4868025 0.4953853 0.5038895 0.5120122 
##    OTU235    OTU121    OTU174     OTU99    OTU925    OTU854    OTU545    OTU144     OTU37 
## 0.5198446 0.5274965 0.5347188 0.5412268 0.5476305 0.5539799 0.5600864 0.5661491 0.5721572 
##    OTU197    OTU140     OTU68     OTU16    OTU340    OTU298    OTU605    OTU640     OTU69 
## 0.5779404 0.5835652 0.5888852 0.5940274 0.5989837 0.6039133 0.6088284 0.6136372 0.6183414 
##    OTU715    OTU105     OTU41    OTU137     OTU88     OTU27   OTU1089     OTU60     OTU45 
## 0.6230415 0.6277280 0.6323451 0.6367965 0.6411466 0.6453214 0.6494391 0.6534993 0.6575361 
##    OTU221    OTU166    OTU218    OTU261    OTU163     OTU81    OTU115    OTU322    OTU259 
## 0.6615384 0.6654843 0.6693905 0.6731251 0.6767907 0.6804246 0.6836906 0.6868900 0.6900614 
##   OTU1047    OTU172    OTU939     OTU67 
## 0.6931820 0.6962866 0.6993799 0.7024051 
## 
## $`Hyperarid loess soil_Hyperarid dolomite`
##      OTU1      OTU2      OTU3      OTU7      OTU6      OTU5      OTU9     OTU18      OTU8 
## 0.1495788 0.2321174 0.3117152 0.3744107 0.4214125 0.4537922 0.4782802 0.5014206 0.5233406 
##     OTU17     OTU41     OTU16    OTU936     OTU11     OTU28     OTU12     OTU20     OTU10 
## 0.5440341 0.5638310 0.5798466 0.5953299 0.6102815 0.6246998 0.6387122 0.6525475 0.6660844 
##     OTU15     OTU14     OTU13     OTU45 
## 0.6762133 0.6855007 0.6942967 0.7030461 
## 
## $`Hyperarid loess soil_Arid loess soil`
##      OTU6     OTU25     OTU62     OTU88     OTU90     OTU11     OTU68    OTU177     OTU78 
## 0.1171964 0.1409598 0.1605516 0.1798885 0.1922238 0.2043007 0.2156943 0.2264171 0.2370192 
##    OTU133    OTU156    OTU187     OTU46    OTU116     OTU67    OTU137    OTU100    OTU687 
## 0.2468319 0.2543986 0.2615122 0.2685441 0.2751451 0.2813878 0.2874536 0.2934299 0.2993293 
##    OTU388    OTU115    OTU121     OTU73    OTU412     OTU91    OTU167     OTU93     OTU16 
## 0.3051050 0.3105215 0.3159249 0.3213181 0.3265045 0.3316752 0.3368401 0.3418967 0.3469359 
##     OTU99    OTU160    OTU218    OTU144    OTU131    OTU135     OTU53    OTU556    OTU152 
## 0.3517643 0.3565360 0.3611535 0.3655646 0.3697266 0.3738519 0.3776564 0.3814527 0.3852354 
##    OTU130    OTU125    OTU227    OTU340    OTU214    OTU233    OTU640    OTU341    OTU311 
## 0.3888725 0.3924354 0.3958973 0.3993476 0.4027687 0.4061872 0.4095931 0.4129866 0.4163480 
##    OTU648    OTU422    OTU716    OTU105    OTU197    OTU163    OTU140     OTU44     OTU37 
## 0.4196998 0.4230194 0.4263048 0.4295726 0.4328230 0.4359927 0.4391577 0.4422934 0.4454224 
##     OTU92    OTU240    OTU715    OTU139    OTU253    OTU221    OTU298    OTU101    OTU171 
## 0.4484404 0.4514561 0.4544639 0.4573782 0.4602606 0.4631238 0.4659602 0.4687870 0.4716058 
##   OTU1235    OTU107    OTU365    OTU491    OTU103    OTU647    OTU479   OTU1047   OTU1288 
## 0.4744165 0.4771309 0.4797889 0.4824360 0.4850823 0.4876817 0.4902569 0.4928293 0.4953810 
##    OTU278    OTU322    OTU136   OTU1102   OTU1050    OTU206    OTU360    OTU220    OTU248 
## 0.4979283 0.5004747 0.5030146 0.5055356 0.5080079 0.5104620 0.5128571 0.5152382 0.5176119 
##    OTU231     OTU48    OTU638    OTU945    OTU315    OTU294    OTU181    OTU190    OTU471 
## 0.5199597 0.5222842 0.5245825 0.5268765 0.5291585 0.5314373 0.5337155 0.5359310 0.5381211 
##    OTU172    OTU169    OTU457     OTU82    OTU191   OTU1013    OTU472    OTU182    OTU251 
## 0.5403020 0.5424802 0.5446509 0.5468014 0.5489293 0.5510565 0.5531619 0.5552609 0.5573352 
##    OTU619    OTU184    OTU352    OTU938    OTU154    OTU707    OTU304    OTU538    OTU386 
## 0.5594017 0.5614639 0.5635207 0.5655628 0.5675880 0.5696106 0.5716293 0.5736343 0.5756356 
##   OTU1116    OTU582    OTU751    OTU434     OTU20   OTU1319   OTU1218    OTU992    OTU309 
## 0.5776195 0.5795812 0.5815142 0.5834436 0.5853509 0.5872476 0.5891360 0.5910228 0.5929026 
##    OTU456    OTU364    OTU302     OTU69    OTU468    OTU725     OTU86    OTU466    OTU562 
## 0.5947562 0.5965751 0.5983597 0.6001340 0.6019002 0.6036514 0.6053964 0.6071295 0.6088489 
##   OTU1184    OTU347    OTU430    OTU164     OTU60    OTU451    OTU635    OTU489    OTU854 
## 0.6105629 0.6122639 0.6139634 0.6156620 0.6173343 0.6189824 0.6206193 0.6222384 0.6238543 
##    OTU610    OTU132    OTU439    OTU255   OTU1046    OTU421    OTU170    OTU543    OTU579 
## 0.6254567 0.6270446 0.6286319 0.6302162 0.6317761 0.6333278 0.6348784 0.6363976 0.6379151 
##   OTU1024     OTU27    OTU323    OTU929    OTU368    OTU953    OTU799    OTU545    OTU258 
## 0.6394210 0.6409169 0.6423915 0.6438630 0.6453326 0.6467955 0.6482429 0.6496722 0.6510934 
##   OTU1282    OTU267     OTU28   OTU1100    OTU277    OTU279    OTU667    OTU199    OTU705 
## 0.6524991 0.6538874 0.6552578 0.6566268 0.6579884 0.6593460 0.6606965 0.6620386 0.6633559 
##    OTU999    OTU232     OTU95    OTU303    OTU836   OTU1242    OTU653      OTU3    OTU141 
## 0.6646697 0.6659821 0.6672939 0.6685879 0.6698524 0.6711151 0.6723738 0.6736310 0.6748814 
##    OTU697    OTU273    OTU499     OTU29   OTU1011    OTU389    OTU162   OTU1108    OTU166 
## 0.6761230 0.6773607 0.6785955 0.6798301 0.6810609 0.6822766 0.6834907 0.6847033 0.6859052 
##    OTU595    OTU524    OTU185    OTU765     OTU41    OTU404    OTU641    OTU557   OTU1005 
## 0.6871070 0.6883061 0.6894992 0.6906885 0.6918492 0.6930093 0.6941662 0.6953181 0.6964525 
##   OTU1049     OTU10   OTU1101    OTU449 
## 0.6975821 0.6987054 0.6998284 0.7009432 
## 
## $`Hyperarid dust_Arid limestone`
##      OTU6     OTU20     OTU40     OTU55    OTU225    OTU422     OTU73     OTU11     OTU46 
## 0.3830030 0.4418277 0.4965700 0.5241168 0.5387383 0.5528407 0.5659891 0.5778580 0.5889906 
##     OTU94    OTU388     OTU29    OTU119     OTU18     OTU34    OTU235     OTU16    OTU854 
## 0.5996827 0.6092874 0.6183609 0.6273489 0.6359934 0.6440744 0.6514858 0.6587135 0.6657503 
##    OTU174    OTU299     OTU68    OTU545      OTU9     OTU37 
## 0.6726643 0.6795149 0.6852951 0.6910019 0.6966084 0.7021179 
## 
## $`Hyperarid dust_Hyperarid dolomite`
##      OTU6      OTU1      OTU2      OTU3      OTU7      OTU5     OTU18      OTU9      OTU8 
## 0.1443995 0.2811765 0.3569750 0.4302446 0.4872402 0.5170148 0.5383456 0.5591551 0.5792831 
##     OTU17     OTU41     OTU11     OTU16    OTU936     OTU28     OTU20     OTU12 
## 0.5983333 0.6167381 0.6326139 0.6480207 0.6622816 0.6756378 0.6885866 0.7014489 
## 
## $`Hyperarid dust_Arid loess soil`
##      OTU6     OTU55     OTU25     OTU88     OTU11     OTU62     OTU68    OTU177     OTU78 
## 0.2717856 0.2922788 0.3102694 0.3274254 0.3433223 0.3590820 0.3688505 0.3770737 0.3852909 
##     OTU94    OTU133    OTU144     OTU46    OTU115     OTU67    OTU156     OTU16    OTU388 
## 0.3931965 0.4006953 0.4078863 0.4149436 0.4216680 0.4280641 0.4338043 0.4395276 0.4449791 
##    OTU116    OTU137    OTU197    OTU100    OTU687     OTU91     OTU73    OTU412    OTU218 
## 0.4500254 0.4550249 0.4599696 0.4648655 0.4695628 0.4742035 0.4785067 0.4824756 0.4864272 
##      OTU9    OTU422    OTU160     OTU53    OTU152    OTU125    OTU135    OTU131    OTU311 
## 0.4902632 0.4938776 0.4974534 0.5009364 0.5043806 0.5077056 0.5109572 0.5141610 0.5171316 
##    OTU556    OTU163    OTU130    OTU484    OTU214    OTU227    OTU240     OTU37    OTU341 
## 0.5200525 0.5228985 0.5256494 0.5283896 0.5310817 0.5337605 0.5364127 0.5390503 0.5416516 
##    OTU233    OTU854    OTU101    OTU648    OTU140    OTU716     OTU92    OTU220    OTU139 
## 0.5442391 0.5467798 0.5493151 0.5518255 0.5543313 0.5568296 0.5592989 0.5617475 0.5641704 
##     OTU44    OTU491     OTU20    OTU278    OTU103    OTU253     OTU60   OTU1235    OTU206 
## 0.5665825 0.5689943 0.5713301 0.5736308 0.5759258 0.5781072 0.5802428 0.5823693 0.5844926 
##    OTU171    OTU107    OTU365    OTU294    OTU647    OTU136    OTU231   OTU1050   OTU1102 
## 0.5865914 0.5886821 0.5907679 0.5928210 0.5948466 0.5968391 0.5988303 0.6007964 0.6027482 
##    OTU479   OTU1288     OTU82    OTU248    OTU434    OTU638    OTU360    OTU945    OTU181 
## 0.6046967 0.6066266 0.6085347 0.6104021 0.6122693 0.6140896 0.6158867 0.6176719 0.6194314 
##    OTU457    OTU472    OTU938   OTU1184    OTU184    OTU315     OTU33    OTU667    OTU471 
## 0.6211887 0.6229334 0.6246741 0.6264112 0.6281455 0.6298719 0.6315934 0.6333143 0.6350343 
##    OTU190   OTU1013    OTU182    OTU538    OTU352    OTU191     OTU86    OTU169    OTU586 
## 0.6367102 0.6383722 0.6400307 0.6416807 0.6433256 0.6449655 0.6465947 0.6482236 0.6498470 
##    OTU304    OTU251    OTU164    OTU619    OTU170    OTU751   OTU1116    OTU707    OTU582 
## 0.6514406 0.6530247 0.6545987 0.6561620 0.6577205 0.6592749 0.6608170 0.6623473 0.6638686 
##    OTU302    OTU154     OTU69    OTU466    OTU725    OTU430    OTU132    OTU545    OTU489 
## 0.6653819 0.6668883 0.6683903 0.6698868 0.6713743 0.6728597 0.6743389 0.6758092 0.6772674 
##    OTU166     OTU48    OTU992    OTU421   OTU1218   OTU1319    OTU456    OTU309    OTU364 
## 0.6787234 0.6801775 0.6816228 0.6830672 0.6845111 0.6859461 0.6873638 0.6887794 0.6901816 
##    OTU635    OTU437    OTU610    OTU562    OTU298    OTU468    OTU347    OTU333 
## 0.6915755 0.6929580 0.6943383 0.6957096 0.6970803 0.6984166 0.6997406 0.7010381 
## 
## $`Arid limestone_Hyperarid dolomite`
##      OTU1      OTU2      OTU3      OTU6      OTU7      OTU5      OTU9      OTU8     OTU18 
## 0.1440378 0.2242599 0.3019930 0.3607941 0.4191722 0.4507340 0.4744385 0.4957463 0.5163290 
##     OTU17     OTU41     OTU40     OTU11    OTU936     OTU16     OTU28     OTU12     OTU20 
## 0.5365375 0.5544474 0.5699925 0.5851556 0.6002748 0.6145844 0.6287286 0.6422857 0.6555544 
##     OTU10     OTU15     OTU14     OTU13     OTU33 
## 0.6687285 0.6786126 0.6876539 0.6962099 0.7045072 
## 
## $`Arid limestone_Arid loess soil`
##       OTU6      OTU40      OTU20      OTU25      OTU88      OTU62      OTU11     OTU225 
## 0.09760135 0.13750741 0.17679560 0.19728621 0.21638051 0.23340169 0.24845894 0.25799747 
##     OTU177      OTU78     OTU133      OTU73      OTU68     OTU119     OTU144      OTU67 
## 0.26735411 0.27646905 0.28479143 0.29307840 0.30008547 0.30708372 0.31394618 0.32068762 
##     OTU156     OTU422      OTU18      OTU29      OTU34     OTU100     OTU116      OTU91 
## 0.32725867 0.33341009 0.33935962 0.34499496 0.35056556 0.35599894 0.36127683 0.36643832 
##     OTU687      OTU46     OTU299     OTU388     OTU235     OTU197     OTU115     OTU218 
## 0.37138684 0.37623323 0.38103236 0.38580627 0.39057213 0.39516382 0.39952480 0.40388431 
##     OTU412     OTU160     OTU174      OTU53     OTU152     OTU925     OTU135     OTU137 
## 0.40821727 0.41229126 0.41633319 0.42024901 0.42416096 0.42806201 0.43172259 0.43528646 
##      OTU37     OTU854     OTU545     OTU163     OTU125     OTU130     OTU131     OTU556 
## 0.43864358 0.44196239 0.44519165 0.44838886 0.45155916 0.45469790 0.45781670 0.46088453 
##     OTU605     OTU311      OTU69     OTU140     OTU341     OTU214      OTU27     OTU233 
## 0.46394974 0.46699883 0.46996880 0.47291024 0.47584516 0.47876919 0.48167258 0.48457197 
##     OTU648     OTU240      OTU16      OTU92      OTU44     OTU491     OTU139     OTU220 
## 0.48743615 0.49029565 0.49310940 0.49591666 0.49871621 0.50143167 0.50413985 0.50682390 
##     OTU227     OTU278     OTU101    OTU1089     OTU103     OTU171      OTU45     OTU253 
## 0.50944639 0.51205855 0.51465647 0.51723675 0.51980481 0.52223855 0.52464256 0.52704332 
##     OTU107     OTU261     OTU365     OTU206     OTU231    OTU1235      OTU41     OTU479 
## 0.52941636 0.53174517 0.53406670 0.53634794 0.53861594 0.54087696 0.54311887 0.54534212 
##    OTU1288    OTU1102      OTU81    OTU1050     OTU294     OTU181      OTU60     OTU945 
## 0.54753509 0.54972106 0.55189155 0.55399097 0.55607199 0.55813344 0.56016996 0.56216653 
##      OTU82      OTU48     OTU360     OTU472     OTU471     OTU136     OTU169     OTU248 
## 0.56416015 0.56614595 0.56813109 0.57010003 0.57205708 0.57400361 0.57593588 0.57786729 
##     OTU315     OTU184     OTU638     OTU182     OTU352     OTU716     OTU647     OTU191 
## 0.57978907 0.58169843 0.58360105 0.58549541 0.58737150 0.58924686 0.59110929 0.59295132 
##     OTU349     OTU457     OTU938    OTU1013     OTU751     OTU619     OTU164     OTU166 
## 0.59479300 0.59662297 0.59843435 0.60023859 0.60201582 0.60378025 0.60554228 0.60727629 
##     OTU251     OTU190     OTU154     OTU434     OTU707    OTU1116    OTU1218     OTU302 
## 0.60900472 0.61070425 0.61239395 0.61407383 0.61574542 0.61739401 0.61904009 0.62068001 
##     OTU489     OTU421     OTU582     OTU286       OTU7     OTU259     OTU667     OTU364 
## 0.62231663 0.62395222 0.62557826 0.62718663 0.62878593 0.63037948 0.63196953 0.63355404 
##     OTU456    OTU1319     OTU298     OTU992     OTU562     OTU304     OTU468     OTU610 
## 0.63512882 0.63670199 0.63825851 0.63979366 0.64131334 0.64282784 0.64433016 0.64582281 
##     OTU170     OTU347     OTU939      OTU86     OTU430     OTU246     OTU239     OTU949 
## 0.64730766 0.64879001 0.65025643 0.65166808 0.65307432 0.65445899 0.65584235 0.65720240 
##     OTU255     OTU439     OTU466     OTU132     OTU530     OTU579     OTU368     OTU362 
## 0.65856178 0.65991755 0.66125549 0.66259309 0.66392616 0.66524781 0.66654856 0.66784058 
##      OTU28     OTU451     OTU653    OTU1299     OTU309     OTU333     OTU399     OTU323 
## 0.66913171 0.67042259 0.67170906 0.67298919 0.67426714 0.67554482 0.67680684 0.67806486 
##     OTU953     OTU348     OTU162    OTU1024     OTU389    OTU1005     OTU929     OTU524 
## 0.67931748 0.68056959 0.68181843 0.68303243 0.68424157 0.68543316 0.68662133 0.68779205 
##    OTU1100     OTU586      OTU95     OTU279     OTU799      OTU87    OTU1101     OTU232 
## 0.68895982 0.69012713 0.69129203 0.69244608 0.69359440 0.69473764 0.69587999 0.69701440 
##     OTU725     OTU303     OTU999 
## 0.69814422 0.69925971 0.70036244 
## 
## $`Hyperarid dolomite_Arid loess soil`
##      OTU1      OTU2      OTU3      OTU7      OTU5      OTU9     OTU18      OTU8     OTU17 
## 0.1319013 0.2062641 0.2786385 0.3307085 0.3600777 0.3818384 0.4015745 0.4212911 0.4401489 
##     OTU41    OTU936     OTU28     OTU12     OTU16     OTU10     OTU11     OTU20      OTU6 
## 0.4571149 0.4711765 0.4842321 0.4967527 0.5091024 0.5210569 0.5325276 0.5434614 0.5531014 
##     OTU15     OTU14     OTU13     OTU24     OTU33     OTU45     OTU88     OTU25     OTU54 
## 0.5622927 0.5706666 0.5785943 0.5859354 0.5932440 0.6005279 0.6077454 0.6147879 0.6216491 
##     OTU62     OTU34     OTU39     OTU47     OTU19     OTU22     OTU35     OTU68     OTU32 
## 0.6280701 0.6340900 0.6399058 0.6455094 0.6510462 0.6560153 0.6601594 0.6640414 0.6677063 
##     OTU29    OTU756     OTU23     OTU37     OTU78     OTU26    OTU177     OTU30    OTU133 
## 0.6713117 0.6748203 0.6783113 0.6817969 0.6852336 0.6886605 0.6920868 0.6953067 0.6984351 
##     OTU46 
## 0.7015026
```

```r
GMPR_ord <- ordinate(Rock_weathering_filt3_GMPR, "CAP", "horn", formula = Rock_weathering_filt3_GMPR ~ Climate * Source)
explained <- eigenvals(GMPR_ord)/sum( eigenvals(GMPR_ord)) * 100
explained <- as.numeric(format(round(explained, 1), nsmall = 1))
data2plot <- cbind(scores(GMPR_ord, display = "sites"), sample_data(Rock_weathering_filt3_GMPR))
p_ord <- ggplot(data2plot) +
  geom_point(aes(x = CAP1, y = CAP2, colour = Source, shape = Climate), size = 3, alpha = 2/3 ) +
  scale_colour_manual(values = pom4) +
  stat_ellipse(aes(x = CAP1, y = CAP2, group = Climate), color = "black", alpha = 0.5, type = "norm", level = 0.95, linetype = 2) +
  xlab(label = paste0("CAP1", " (", explained[1],"%)")) + 
  ylab(label = paste0("CAP2", " (", explained[2],"%)")) +
  coord_fixed(sqrt(explained[2] / explained[1])) 
  
print(p_ord)
```

![](Rock_weathering_figures/ordination final-1.svg)<!-- -->

### Taxonomic features
Explore and plot the taxonomic distribution of the sequences


```r
Rock_weathering_filt3_GMPR_rel <- transform_sample_counts(Rock_weathering_filt3_GMPR, function(x) x / sum(x) )

Rock_weathering_filt3_GMPR_rel %>% 
  sample_data() %>% 
  arrange(Climate, Source) %>% 
  .$sample_names ->
  Sample_order

Rock_weathering_filt3_100 <-
  prune_taxa(names(sort(taxa_sums(Rock_weathering_filt3_GMPR_rel), TRUE)[1:100]), Rock_weathering_filt3_GMPR_rel)
plot_heatmap(
  Rock_weathering_filt3_100,
  method = NULL,
  distance = NULL,
  sample.label = "sample_title",
  taxa.label = "Order",
  sample_order = Sample_order,
  low = "#000033",
  high = "#FF3300"
) #+ theme_bw(base_size = 20) + theme(axis.text.x = element_text(hjust = 0, angle = -90.0))
```

![](Rock_weathering_figures/seqs heatmaps-1.svg)<!-- -->

Let's look at the agglomerated taxa


```r
Rock_weathering_filt3_glom <- tax_glom(Rock_weathering_filt3_GMPR, 
                             "Phylum", 
                             NArm = TRUE)
Rock_weathering_filt3_glom_rel <- transform_sample_counts(Rock_weathering_filt3_glom, function(x) x / sum(x)) 
Rock_weathering_filt3_glom_rel_DF <- psmelt(Rock_weathering_filt3_glom_rel)
Rock_weathering_filt3_glom_rel_DF$Phylum %<>% as.character()

# group dataframe by Phylum, calculate median rel. abundance
Rock_weathering_filt3_glom_rel_DF %>%
  group_by(Phylum) %>%
  summarise(median = median(Abundance)) ->
  medians

# find Phyla whose rel. abund. is less than 0.5%
Rare_phyla <- medians[medians$median <= 0.005, ]$Phylum

# change their name to "Rare"
Rock_weathering_filt3_glom_rel_DF[Rock_weathering_filt3_glom_rel_DF$Phylum %in% Rare_phyla, ]$Phylum <- 'Rare'
# re-group
Rock_weathering_filt3_glom_rel_DF %>%
  group_by(Sample, Climate, Phylum, Rock.type, Source) %>%
  summarise(Abundance = sum(Abundance)) ->
  Rock_weathering_filt3_glom_rel_DF_2plot

# ab.taxonomy$Freq <- sqrt(ab.taxonomy$Freq)
Rock_weathering_filt3_glom_rel_DF_2plot$Phylum %<>% sub("unclassified", "Unclassified", .)
Rock_weathering_filt3_glom_rel_DF_2plot$Phylum %<>% sub("uncultured", "Unclassified", .)

Rock_weathering_filt3_glom_rel_DF_2plot %>% 
  group_by(Sample) %>% 
  filter(Phylum == "Rare") %>% 
  summarise(`Rares (%)` = sum(Abundance * 100)) -> 
  Rares
# Percentage of reads classified as rare 
Rares %>% 
  kable(., digits = 2, caption = "Percentage of reads per sample type classified as rare:") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Percentage of reads per sample type classified as rare:</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:right;"> Rares (%) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SbDust1S14 </td>
   <td style="text-align:right;"> 2.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbDust2S31 </td>
   <td style="text-align:right;"> 0.43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp1SNW49 </td>
   <td style="text-align:right;"> 0.82 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp2SNW50 </td>
   <td style="text-align:right;"> 0.91 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp3SNW51 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp4SNW52 </td>
   <td style="text-align:right;"> 1.99 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp5SNW53 </td>
   <td style="text-align:right;"> 0.72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp6SNW54 </td>
   <td style="text-align:right;"> 1.59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSoil1SA10 </td>
   <td style="text-align:right;"> 9.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSoil2SA11 </td>
   <td style="text-align:right;"> 10.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSoil3SA12 </td>
   <td style="text-align:right;"> 9.62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad1SNW55 </td>
   <td style="text-align:right;"> 2.90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad2SNW56 </td>
   <td style="text-align:right;"> 0.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad3SNW57 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad4SNW58 </td>
   <td style="text-align:right;"> 0.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad5SNW59 </td>
   <td style="text-align:right;"> 0.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad6SNW60 </td>
   <td style="text-align:right;"> 1.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvDust1S32 </td>
   <td style="text-align:right;"> 1.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvDust2S33 </td>
   <td style="text-align:right;"> 0.12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp1GS70 </td>
   <td style="text-align:right;"> 1.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp2GS71 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp3CS25 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp3GS72 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp4GS73 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp5GS74 </td>
   <td style="text-align:right;"> 0.72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp6GS75 </td>
   <td style="text-align:right;"> 0.34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad1GS76 </td>
   <td style="text-align:right;"> 9.37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad2CS23 </td>
   <td style="text-align:right;"> 2.17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad2GS77 </td>
   <td style="text-align:right;"> 0.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad3CS27 </td>
   <td style="text-align:right;"> 2.96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad3GS78 </td>
   <td style="text-align:right;"> 0.16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad4GS79 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad5GS80 </td>
   <td style="text-align:right;"> 0.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad6GS81 </td>
   <td style="text-align:right;"> 0.38 </td>
  </tr>
</tbody>
</table>

```r
sample_order <- match(Rares$Sample, row.names(sample_data(Rock_weathering_filt3_glom)))
Rares %<>% arrange(., sample_order)

Rares %>% 
  cbind(., sample_data(Rock_weathering_filt3_glom)) %>% 
  group_by(Climate.Source) %>% 
  summarise(`Rares (%)` = mean(`Rares (%)`)) -> 
  Rares_merged

# Percentage of reads classified as rare 
Rares %>% 
  kable(., digits = 2, caption = "Percentage of reads per sample type classified as rare:") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Percentage of reads per sample type classified as rare:</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:right;"> Rares (%) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SbDust1S14 </td>
   <td style="text-align:right;"> 2.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad2CS23 </td>
   <td style="text-align:right;"> 2.17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp3CS25 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad3CS27 </td>
   <td style="text-align:right;"> 2.96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbDust2S31 </td>
   <td style="text-align:right;"> 0.43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvDust1S32 </td>
   <td style="text-align:right;"> 1.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvDust2S33 </td>
   <td style="text-align:right;"> 0.12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp1SNW49 </td>
   <td style="text-align:right;"> 0.82 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp2SNW50 </td>
   <td style="text-align:right;"> 0.91 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp3SNW51 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp4SNW52 </td>
   <td style="text-align:right;"> 1.99 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp5SNW53 </td>
   <td style="text-align:right;"> 0.72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSlp6SNW54 </td>
   <td style="text-align:right;"> 1.59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad1SNW55 </td>
   <td style="text-align:right;"> 2.90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad2SNW56 </td>
   <td style="text-align:right;"> 0.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad3SNW57 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad4SNW58 </td>
   <td style="text-align:right;"> 0.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad5SNW59 </td>
   <td style="text-align:right;"> 0.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbWad6SNW60 </td>
   <td style="text-align:right;"> 1.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp1GS70 </td>
   <td style="text-align:right;"> 1.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp2GS71 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp3GS72 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp4GS73 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp5GS74 </td>
   <td style="text-align:right;"> 0.72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvSlp6GS75 </td>
   <td style="text-align:right;"> 0.34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad1GS76 </td>
   <td style="text-align:right;"> 9.37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad2GS77 </td>
   <td style="text-align:right;"> 0.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad3GS78 </td>
   <td style="text-align:right;"> 0.16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad4GS79 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad5GS80 </td>
   <td style="text-align:right;"> 0.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UvWad6GS81 </td>
   <td style="text-align:right;"> 0.38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSoil1SA10 </td>
   <td style="text-align:right;"> 9.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSoil2SA11 </td>
   <td style="text-align:right;"> 10.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SbSoil3SA12 </td>
   <td style="text-align:right;"> 9.62 </td>
  </tr>
</tbody>
</table>

```r
Rock_weathering_filt3_glom_rel_DF_2plot %>% 
  group_by(Phylum) %>% 
  summarise(sum.Taxa = sum(Abundance)) %>% 
  arrange(desc(sum.Taxa)) -> Taxa_rank
Rock_weathering_filt3_glom_rel_DF_2plot$Phylum %<>% 
  factor(., levels = Taxa_rank$Phylum) %>% 
  fct_relevel(., "Rare", after = Inf)
  
p_taxa_box <-
  ggplot(Rock_weathering_filt3_glom_rel_DF_2plot, aes(x = Phylum, y = (Abundance * 100))) +
  geom_boxplot(aes(group = interaction(Phylum, Source)), position = position_dodge(width = 0.9), fatten = 1) +
  geom_point(
    aes(colour = Source),
    position = position_jitterdodge(dodge.width = 1),
    alpha = 1 / 2,
    stroke = 0,
    size = 2
  ) +
  scale_colour_manual(values = pom4, name = "") +
  theme_cowplot(font_size = 11, font_family = f_name) +
  labs(x = NULL, y = "Relative abundance (%)") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  facet_grid(Climate ~ .) +
  background_grid(major = "xy",
                  minor = "none") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.9,
    hjust = 0.9
  ))
print(p_taxa_box)
```

![](Rock_weathering_figures/agglomerated taxa box plot-1.svg)<!-- -->

#### Test differences between samples on the phylum level

```r
Taxa_tests_phylum <- STAMPR(Rock_weathering_filt3_GMPR, "Phylum", sig_pairs)
```

```
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Phylum          
## OTU1 "Proteobacteria"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1 1118.38 11.2778 0.00078
## Source          3  412.13  4.1560 0.24511
## Climate:Source  1    2.40  0.0242 0.87637
## Residuals      28 1739.58                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -57.31870  77.10514
## sample estimates:
## difference in location 
##              -8.837742 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 13, p-value = 0.9273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -45.05308  61.68040
## sample estimates:
## difference in location 
##               3.711094 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 3, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -44.64033 -10.50860
## sample estimates:
## difference in location 
##              -24.33348 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -31.61913  23.25576
## sample estimates:
## difference in location 
##               2.734424 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -92.69994 -22.88429
## sample estimates:
## difference in location 
##               -41.7985 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.006099
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -43.71818 -14.50227
## sample estimates:
## difference in location 
##              -31.94182 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 12, p-value = 0.4273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -9.149286  9.089127
## sample estimates:
## difference in location 
##              -3.911288 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 3, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -43.181793  -5.626402
## sample estimates:
## difference in location 
##              -35.79785 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Phylum               
## OTU2 "Deinococcus-Thermus"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df Sum Sq       H p.value
## Climate         1  212.5  2.1429 0.14323
## Source          3 1474.7 14.8713 0.00193
## Climate:Source  1   72.6  0.7321 0.39220
## Residuals      28 1512.7                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.1159525 28.8971507
## sample estimates:
## difference in location 
##               3.763569 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -16.468193  -6.626912
## sample estimates:
## difference in location 
##              -12.42899 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.035956 20.819838
## sample estimates:
## difference in location 
##              0.1319845 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 31, p-value = 0.0712
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.559324 12.947540
## sample estimates:
## difference in location 
##               6.279289 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   6.553372 16.590053
## sample estimates:
## difference in location 
##               11.95774 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 53, p-value = 0.2855
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -11.37397   7.94874
## sample estimates:
## difference in location 
##              -5.992007 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 35, p-value = 0.01724
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.1726432 24.9922558
## sample estimates:
## difference in location 
##               2.545225 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -16.14293  -6.64311
## sample estimates:
## difference in location 
##              -11.97071 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Phylum         
## OTU5 "Bacteroidetes"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H p.value
## Climate         1   12.97 0.1308 0.71761
## Source          3  821.88 8.2879 0.04042
## Climate:Source  1   74.82 0.7545 0.38507
## Residuals      28 2362.83               
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -68.34270  13.97676
## sample estimates:
## difference in location 
##              -27.18296 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 16, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -10.44161  67.94132
## sample estimates:
## difference in location 
##               29.64229 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 26, p-value = 0.279
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.354092  4.988486
## sample estimates:
## difference in location 
##               2.066666 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 28, p-value = 0.1703
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.9948177 10.0975854
## sample estimates:
## difference in location 
##                5.62787 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.8221862 15.0872946
## sample estimates:
## difference in location 
##               7.651679 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 51, p-value = 0.2366
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -6.951583  1.649500
## sample estimates:
## difference in location 
##              -3.570295 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 10, p-value = 0.279
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -5.569969  2.833986
## sample estimates:
## difference in location 
##              -1.641379 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 15, p-value = 0.7182
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.943049  5.141081
## sample estimates:
## difference in location 
##              -2.618944 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Phylum         
## OTU7 "Cyanobacteria"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H p.value
## Climate         1   95.56 0.9636 0.32628
## Source          3  853.37 8.6055 0.03502
## Climate:Source  1   45.07 0.4545 0.50023
## Residuals      28 2278.50               
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 19, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.455407 18.180308
## sample estimates:
## difference in location 
##               2.957942 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 6, p-value = 0.3153
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -37.738611   2.344579
## sample estimates:
## difference in location 
##              -3.324833 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 33, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.6713462 18.0874183
## sample estimates:
## difference in location 
##               3.080276 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 34, p-value = 0.02527
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.06861458 24.43779425
## sample estimates:
## difference in location 
##               5.272126 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 20, p-value = 0.1709
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5906178 37.8005207
## sample estimates:
## difference in location 
##               4.208215 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 73, p-value = 0.977
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.311598  5.103880
## sample estimates:
## difference in location 
##              0.2627632 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.7182
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.692326 14.614652
## sample estimates:
## difference in location 
##                1.12547 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 16, p-value = 0.8286
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -23.702686   4.781333
## sample estimates:
## difference in location 
##              -1.638922 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Phylum         
## OTU10 "Acidobacteria"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1    2.38  0.0240 0.87682
## Source          3 1299.20 13.1012 0.00442
## Climate:Source  1   70.42  0.7101 0.39942
## Residuals      28 1900.50                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 10, p-value = 0.7842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.875998  6.843540
## sample estimates:
## difference in location 
##             -0.2000774 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.891498  1.987482
## sample estimates:
## difference in location 
##              -1.607997 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.440454  1.139749
## sample estimates:
## difference in location 
##           -0.007451021 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 28, p-value = 0.1703
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.8000388  2.8098831
## sample estimates:
## difference in location 
##                1.66951 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.1207
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.8738157  4.2049748
## sample estimates:
## difference in location 
##               2.407387 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 40, p-value = 0.06896
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.5382471  0.1114779
## sample estimates:
## difference in location 
##              -1.540419 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -10.418316  -7.394329
## sample estimates:
## difference in location 
##               -8.59669 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  5.710241 9.481740
## sample estimates:
## difference in location 
##               7.038468 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Phylum          
## OTU20 "Actinobacteria"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H p.value
## Climate         1  584.74 5.8965 0.01517
## Source          3  761.45 7.6785 0.05315
## Climate:Source  1    0.82 0.0082 0.92769
## Residuals      28 1925.50               
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -10.12738  56.93452
## sample estimates:
## difference in location 
##               36.90999 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 7, p-value = 0.4113
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -44.932741   9.863263
## sample estimates:
## difference in location 
##              -6.636044 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 29, p-value = 0.1296
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -6.174561 41.738766
## sample estimates:
## difference in location 
##               20.55627 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -23.05721  22.19200
## sample estimates:
## difference in location 
##              -6.063935 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.05523
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.954275 46.692652
## sample estimates:
## difference in location 
##               12.84274 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 115, p-value = 0.01414
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   7.16780 40.31334
## sample estimates:
## difference in location 
##                21.9016 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 27, p-value = 0.2199
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -11.03421  33.22758
## sample estimates:
## difference in location 
##               17.11401 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 25, p-value = 0.3481
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -19.10688  25.75876
## sample estimates:
## difference in location 
##               15.87242 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Phylum      
## OTU44 "Firmicutes"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  362.38  3.6543 0.05593
## Source          3 1569.80 15.8299 0.00123
## Climate:Source  1    1.07  0.0108 0.91740
## Residuals      28 1339.25                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.618072 -4.812207
## sample estimates:
## difference in location 
##              -6.671831 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.2551915 8.7211006
## sample estimates:
## difference in location 
##               5.343429 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 2, p-value = 0.02527
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.01263008 -0.03357744
## sample estimates:
## difference in location 
##             -0.3462271 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.516
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.6665798  2.1292269
## sample estimates:
## difference in location 
##              0.1278182 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 8, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.274515  4.521171
## sample estimates:
## difference in location 
##             -0.2052955 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 20, p-value = 0.002946
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.6614044 -0.2036041
## sample estimates:
## difference in location 
##             -0.5325774 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.2671454 -0.7882357
## sample estimates:
## difference in location 
##             -0.9980437 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.505982  1.119910
## sample estimates:
## difference in location 
##              0.5223647 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##        Phylum            
## OTU116 "Gemmatimonadetes"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  536.03  5.4053 0.02008
## Source          3 1335.07 13.4629 0.00374
## Climate:Source  1   29.40  0.2965 0.58610
## Residuals      28 1372.00                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 13, p-value = 0.9273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.563363  2.368097
## sample estimates:
## difference in location 
##              0.1096169 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.7842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5939111  1.5794065
## sample estimates:
## difference in location 
##              0.4927523 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.0712
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.6300118  0.1465505
## sample estimates:
## difference in location 
##              -1.280392 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.8652168 -0.9657726
## sample estimates:
## difference in location 
##              -1.923801 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.1207
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.1146488  0.6986154
## sample estimates:
## difference in location 
##              0.2262944 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 115, p-value = 0.01414
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.06797809 1.12543652
## sample estimates:
## difference in location 
##               0.575367 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -5.196120 -3.074545
## sample estimates:
## difference in location 
##              -4.256295 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  4.371457 5.332202
## sample estimates:
## difference in location 
##               4.958494 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##        Phylum       
## OTU225 "Chloroflexi"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  724.97  7.3106 0.00685
## Source          3 1423.13 14.3509 0.00246
## Climate:Source  1   17.07  0.1721 0.67825
## Residuals      28 1107.33                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.888273 27.556451
## sample estimates:
## difference in location 
##                8.68754 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 13, p-value = 0.9273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.388231  2.129751
## sample estimates:
## difference in location 
##              0.2617992 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 6, p-value = 0.09694
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -19.842294   5.726156
## sample estimates:
## difference in location 
##              -8.397385 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -23.17677 -10.76493
## sample estimates:
## difference in location 
##              -21.05372 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.7842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.822953  2.765635
## sample estimates:
## difference in location 
##              0.1344354 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 133, p-value = 0.0004777
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   4.219769 15.572571
## sample estimates:
## difference in location 
##               8.960243 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 10, p-value = 0.279
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -15.090981   6.626764
## sample estimates:
## difference in location 
##              -6.713824 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  13.68642 17.40189
## sample estimates:
## difference in location 
##                16.2042
```

```r
Taxa_tests_order <- STAMPR(Rock_weathering_filt3_GMPR, "Order", sig_pairs)
```

```
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Order            
## OTU1 "Burkholderiales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  900.74  9.0830 0.00258
## Source          3 1292.08 13.0294 0.00457
## Climate:Source  1   79.35  0.8002 0.37104
## Residuals      28 1000.33                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -41.194215  -3.356289
## sample estimates:
## difference in location 
##              -22.27526 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 12, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -40.18183  41.15491
## sample estimates:
## difference in location 
##             -0.3563909 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 2, p-value = 0.02527
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.4548998 -0.1261025
## sample estimates:
## difference in location 
##             -0.2969459 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 33, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  12.65335 40.77171
## sample estimates:
## difference in location 
##               34.15781 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.325322 43.909135
## sample estimates:
## difference in location 
##                34.7954 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 9, p-value = 0.000308
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -36.54463 -15.27047
## sample estimates:
## difference in location 
##              -34.67256 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.1343499 -0.3934692
## sample estimates:
## difference in location 
##               -0.74279 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 3, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -40.52361 -11.97382
## sample estimates:
## difference in location 
##              -33.78455 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Order               
## OTU5 "Sphingobacteriales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  183.56  1.8510 0.17367
## Source          3 1245.17 12.5564 0.00570
## Climate:Source  1   11.27  0.1136 0.73607
## Residuals      28 1832.50                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 8, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.1104079  0.8423098
## sample estimates:
## difference in location 
##             -0.3089557 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -13.082829   1.131093
## sample estimates:
## difference in location 
##              -6.123185 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.7182
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.6640484  0.7639180
## sample estimates:
## difference in location 
##             0.01919754 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 30, p-value = 0.09694
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5096519 10.2461215
## sample estimates:
## difference in location 
##               6.491376 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.1207
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.4782018 13.6953205
## sample estimates:
## difference in location 
##                6.46166 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.005108
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.5172328 -0.8528489
## sample estimates:
## difference in location 
##              -6.135539 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.329414 -1.612209
## sample estimates:
## difference in location 
##              -3.002078 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.039024  3.572766
## sample estimates:
## difference in location 
##              -3.718278 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Order              
## OTU6 "Enterobacteriales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1   70.62  0.7121 0.39874
## Source          3 2314.23 23.3368 0.00003
## Climate:Source  1   64.07  0.6461 0.42153
## Residuals      28  823.58                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 20, p-value = 0.1709
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5713626 91.1323680
## sample estimates:
## difference in location 
##               3.553661 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -10.437379   1.912076
## sample estimates:
## difference in location 
##                1.69055 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 3, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -49.501489  -5.472717
## sample estimates:
## difference in location 
##              -21.46237 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -51.68982 -22.42206
## sample estimates:
## difference in location 
##              -26.58309 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -84.18527 -55.63994
## sample estimates:
## difference in location 
##              -69.91261 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 134, p-value = 0.0003842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  1.932528 8.058776
## sample estimates:
## difference in location 
##                5.11499 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   1.298489 17.365156
## sample estimates:
## difference in location 
##               5.334665 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 8, p-value = 0.1703
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.58435850  0.06542455
## sample estimates:
## difference in location 
##             -0.1028499 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Order         
## OTU7 "SubsectionII"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H p.value
## Climate         1  174.38 1.7585 0.18481
## Source          3  716.55 7.2257 0.06504
## Climate:Source  1    2.40 0.0242 0.87637
## Residuals      28 2379.17               
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.7842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.04641 14.86629
## sample estimates:
## difference in location 
##              0.1945796 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 4, p-value = 0.1709
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -33.5566632   0.8541497
## sample estimates:
## difference in location 
##              -3.536299 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.2696644  9.1526642
## sample estimates:
## difference in location 
##              0.5183234 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 34, p-value = 0.02527
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.1576628 23.9114605
## sample estimates:
## difference in location 
##               4.181034 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.1207
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.3651775 33.6630338
## sample estimates:
## difference in location 
##               4.549962 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 41, p-value = 0.07825
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.940779  1.057636
## sample estimates:
## difference in location 
##              -1.222737 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 20, p-value = 0.8286
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5775934  9.0335277
## sample estimates:
## difference in location 
##              0.2127829 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.0712
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -23.7924121   0.1227166
## sample estimates:
## difference in location 
##              -3.891536 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##      Order        
## OTU9 "Rhizobiales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H  p.value
## Climate         1  880.26 8.8766 0.002888
## Source          3  546.82 5.5141 0.137795
## Climate:Source  1  400.42 4.0378 0.044491
## Residuals      28 1445.00                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 9, p-value = 0.6481
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -5.881060  2.947064
## sample estimates:
## difference in location 
##              -1.466998 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 6, p-value = 0.3153
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -11.662705   5.901289
## sample estimates:
## difference in location 
##              -3.357776 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 16, p-value = 0.8286
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.990814  1.211501
## sample estimates:
## difference in location 
##             -0.5011302 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 33, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.9757943 10.1285233
## sample estimates:
## difference in location 
##               4.438603 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 8, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.848487  4.603742
## sample estimates:
## difference in location 
##              -1.200014 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 12, p-value = 0.000592
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.323014 -2.448342
## sample estimates:
## difference in location 
##                -4.6125 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.588206 -2.234018
## sample estimates:
## difference in location 
##              -2.956213 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 12, p-value = 0.4273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.265095  1.276736
## sample estimates:
## difference in location 
##              -1.607365 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order       
## OTU10 "Subgroup_4"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1   49.44  0.4986 0.48013
## Source          3 1314.88 13.2592 0.00411
## Climate:Source  1   72.60  0.7321 0.39220
## Residuals      28 1835.58                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 10, p-value = 0.7842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.355819  6.528750
## sample estimates:
## difference in location 
##             -0.1126584 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 4, p-value = 0.1709
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.075922  1.321435
## sample estimates:
## difference in location 
##              -1.911203 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.9403130  0.8227086
## sample estimates:
## difference in location 
##              -0.265193 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 30, p-value = 0.09694
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5559751  3.2464937
## sample estimates:
## difference in location 
##               1.832547 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.6122869  4.2773329
## sample estimates:
## difference in location 
##               2.400129 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 31, p-value = 0.01937
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.8472840 -0.3238103
## sample estimates:
## difference in location 
##              -1.969475 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 3, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -6.169976 -4.029602
## sample estimates:
## difference in location 
##              -5.295554 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  1.716829 5.439081
## sample estimates:
## difference in location 
##               3.039978 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order             
## OTU11 "Sphingomonadales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H p.value
## Climate         1    4.97 0.0501 0.82285
## Source          3  595.88 6.0089 0.11118
## Climate:Source  1   66.15 0.6671 0.41408
## Residuals      28 2605.50               
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 13, p-value = 0.9273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -5.903441 12.844720
## sample estimates:
## difference in location 
##                0.13997 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 9, p-value = 0.6481
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.326224  4.103880
## sample estimates:
## difference in location 
##             -0.7099819 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 20, p-value = 0.8286
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.128869  5.956528
## sample estimates:
## difference in location 
##              0.3258721 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.6134
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.226862  5.505896
## sample estimates:
## difference in location 
##               1.482761 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  1.378665 9.581825
## sample estimates:
## difference in location 
##               5.117695 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 52, p-value = 0.2602
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.503731  1.094646
## sample estimates:
## difference in location 
##              -1.199836 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 12, p-value = 0.4273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.995816  3.069399
## sample estimates:
## difference in location 
##                -1.0588 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 15, p-value = 0.7182
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.890330  2.811984
## sample estimates:
## difference in location 
##             -0.5360087 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order            
## OTU14 "Caulobacterales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  430.62  4.3424 0.03718
## Source          3 1247.53 12.5802 0.00564
## Climate:Source  1   40.02  0.4035 0.52527
## Residuals      28 1554.33                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 9, p-value = 0.6481
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.0704490  0.9475039
## sample estimates:
## difference in location 
##            -0.03852294 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.9373695  0.5232694
## sample estimates:
## difference in location 
##             -0.8493462 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.4273
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.4313298  0.7733929
## sample estimates:
## difference in location 
##               0.218237 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 35, p-value = 0.01724
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.2813935 2.4611436
## sample estimates:
## difference in location 
##               1.096272 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.3652594 3.0358987
## sample estimates:
## difference in location 
##               1.288226 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 17, p-value = 0.001652
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.6797914 -0.4360266
## sample estimates:
## difference in location 
##             -0.9214106 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.4322340  0.3639032
## sample estimates:
## difference in location 
##             0.01551463 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.043705 -0.096489
## sample estimates:
## difference in location 
##             -0.9607206 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order            
## OTU20 "Rubrobacterales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  362.38  3.6543 0.05593
## Source          3 1283.60 12.9439 0.00476
## Climate:Source  1    7.35  0.0741 0.78543
## Residuals      28 1619.17                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   1.343966 35.280345
## sample estimates:
## difference in location 
##               23.05921 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -34.39409979  -0.07048124
## sample estimates:
## difference in location 
##              -4.199694 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 33, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   5.206833 27.553233
## sample estimates:
## difference in location 
##               17.11954 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 19, p-value = 0.9425
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -6.01319 27.58079
## sample estimates:
## difference in location 
##              0.3061655 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.05523
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.1729386 34.3585874
## sample estimates:
## difference in location 
##               4.639865 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 110, p-value = 0.03038
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   1.007057 21.596004
## sample estimates:
## difference in location 
##               13.50276 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 33, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   5.377442 24.819725
## sample estimates:
## difference in location 
##               16.82808 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 20, p-value = 0.8286
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -26.414743   5.855238
## sample estimates:
## difference in location 
##                1.69516 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order         
## OTU22 "Cytophagales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1 1164.74 11.7452 0.00061
## Source          3   44.85  0.4522 0.92925
## Climate:Source  1  281.67  2.8403 0.09192
## Residuals      28 1781.25                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 7, p-value = 0.4113
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -67.49137  13.79442
## sample estimates:
## difference in location 
##              -27.69087 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.9995852 67.2950902
## sample estimates:
## difference in location 
##               34.90214 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 26, p-value = 0.279
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.069587  5.465177
## sample estimates:
## difference in location 
##               1.337481 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 10, p-value = 0.279
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.8392488  0.3526159
## sample estimates:
## difference in location 
##             -0.2786225 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.05523
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.1948212  1.8870346
## sample estimates:
## difference in location 
##              0.5415788 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 115, p-value = 0.01414
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.416069 3.919566
## sample estimates:
## difference in location 
##               2.374847 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.516
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.866692  5.145857
## sample estimates:
## difference in location 
##               1.171801 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 34, p-value = 0.02527
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.1545652 1.8629482
## sample estimates:
## difference in location 
##               1.157305 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order          
## OTU30 "Micrococcales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq      H p.value
## Climate         1  265.44 2.6767 0.10183
## Source          3  749.89 7.5619 0.05599
## Climate:Source  1   33.75 0.3403 0.55964
## Residuals      28 2223.42               
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 2, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -5.206538  4.628937
## sample estimates:
## difference in location 
##              -2.182672 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.828908  5.114544
## sample estimates:
## difference in location 
##               2.137151 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 25, p-value = 0.3481
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5320493  1.4765078
## sample estimates:
## difference in location 
##              0.2706085 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 26, p-value = 0.279
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5122759  0.7866787
## sample estimates:
## difference in location 
##              0.2027858 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 7, p-value = 0.4113
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.5000857  5.5661994
## sample estimates:
## difference in location 
##             -0.2281722 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 79, p-value = 0.7075
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.2204756  0.3690566
## sample estimates:
## difference in location 
##              0.0729464 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 4, p-value = 0.05135
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.5258004  0.2675506
## sample estimates:
## difference in location 
##              -1.065072 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 33, p-value = 0.03636
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.4223112 1.4948001
## sample estimates:
## difference in location 
##               1.103329 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order       
## OTU33 "Frankiales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1   64.97 0.65517 0.41827
## Source          3  262.68 2.64887 0.44899
## Climate:Source  1    9.60 0.09681 0.75570
## Residuals      28 2935.25                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 16, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.724773 10.232076
## sample estimates:
## difference in location 
##               1.105793 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 8, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -6.780411  2.832335
## sample estimates:
## difference in location 
##             -0.4845265 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 19, p-value = 0.9425
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.568094  4.262904
## sample estimates:
## difference in location 
##              0.2980444 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 16, p-value = 0.8286
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.577901  4.274808
## sample estimates:
## difference in location 
##             -0.1588814 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 0.3153
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.7863939  6.5607849
## sample estimates:
## difference in location 
##               1.092154 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 82, p-value = 0.5834
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.373038  1.650491
## sample estimates:
## difference in location 
##              0.5690406 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 17, p-value = 0.9425
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.015699  2.983868
## sample estimates:
## difference in location 
##             -0.1357949 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.7182
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.105227  1.803981
## sample estimates:
## difference in location 
##              0.4714685 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order          
## OTU40 "Deinococcales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1   24.74  0.2494 0.61748
## Source          3 1445.20 14.5734 0.00222
## Climate:Source  1  190.82  1.9242 0.16539
## Residuals      28 1611.75                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.1159525 28.8971507
## sample estimates:
## difference in location 
##               3.763569 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.3472144 -0.5879807
## sample estimates:
## difference in location 
##              -2.500071 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.04070 20.81984
## sample estimates:
## difference in location 
##              0.1319845 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 11, p-value = 0.3481
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -7.129548  2.026183
## sample estimates:
## difference in location 
##              -2.583399 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 24, p-value = 0.03576
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.5144736 7.4691341
## sample estimates:
## difference in location 
##               2.253139 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 82, p-value = 0.5834
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -1.344459 18.286424
## sample estimates:
## difference in location 
##              0.9522018 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   0.3237722 25.0161232
## sample estimates:
## difference in location 
##               2.685594 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.2047144 -0.6847059
## sample estimates:
## difference in location 
##              -2.525927 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##       Order                
## OTU73 "Solirubrobacterales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  670.62  6.7625 0.00931
## Source          3 1220.47 12.3072 0.00640
## Climate:Source  1    3.75  0.0378 0.84581
## Residuals      28 1377.67                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.2353533 27.6409333
## sample estimates:
## difference in location 
##               4.612546 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 5, p-value = 0.2353
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -2.9351537  0.2608993
## sample estimates:
## difference in location 
##             -0.4958395 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.212403  9.958079
## sample estimates:
## difference in location 
##             -0.1090627 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -5.293649 -2.537352
## sample estimates:
## difference in location 
##              -4.224644 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 21, p-value = 0.1207
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.1764773  3.0664060
## sample estimates:
## difference in location 
##              0.5754338 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 125, p-value = 0.002437
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##   1.196491 11.780432
## sample estimates:
## difference in location 
##               3.693183 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 17, p-value = 0.9425
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.307577 10.065980
## sample estimates:
## difference in location 
##             -0.4198441 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  2.036055 5.953411
## sample estimates:
## difference in location 
##               3.962095 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##        Order             
## OTU144 "Acidimicrobiales"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1  402.62  4.0600 0.04391
## Source          3 1958.32 19.7477 0.00019
## Climate:Source  1   36.82  0.3713 0.54232
## Residuals      28  874.75                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.2770006  1.8753928
## sample estimates:
## difference in location 
##              0.9466836 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 8, p-value = 0.5228
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.6952939  0.2939760
## sample estimates:
## difference in location 
##            -0.07693071 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.055194 -2.733487
## sample estimates:
## difference in location 
##              -4.376145 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -8.811553 -3.722903
## sample estimates:
## difference in location 
##              -5.306655 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 23, p-value = 0.05523
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.009727731  0.707286393
## sample estimates:
## difference in location 
##              0.2041089 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 125, p-value = 0.002437
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  0.3522632 1.1649116
## sample estimates:
## difference in location 
##              0.8707461 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.351919 -1.703828
## sample estimates:
## difference in location 
##              -2.347733 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  2.648734 4.108222
## sample estimates:
## difference in location 
##               3.124562 
## 
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##        Order         
## OTU225 "JG30-KF-CM45"
## 
## DV:  Abundance 
## Observations:  34 
## D:  1 
## MS total:  99.16667 
## 
##                Df  Sum Sq       H p.value
## Climate         1 1050.62 10.5945 0.00113
## Source          3 1025.72 10.3434 0.01586
## Climate:Source  1   20.42  0.2059 0.65001
## Residuals      28 1175.75                
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 22, p-value = 0.08284
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.3350587 20.3127324
## sample estimates:
## difference in location 
##               4.143156 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 14, p-value = 0.7842
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.4073917  0.4038865
## sample estimates:
## difference in location 
##            0.007441521 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 19, p-value = 0.9425
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.843441 11.810661
## sample estimates:
## difference in location 
##              0.6013048 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 0, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -4.795687 -2.220769
## sample estimates:
## difference in location 
##              -4.477859 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 18, p-value = 0.3153
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.1554704  0.5143634
## sample estimates:
## difference in location 
##              0.1156305 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 135, p-value = 0.000308
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  2.579639 6.320156
## sample estimates:
## difference in location 
##               4.241255 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 17, p-value = 0.9425
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -3.807166 10.334104
## sample estimates:
## difference in location 
##             -0.2087371 
## 
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  Abundance by Climate.Source
## W = 36, p-value = 0.01154
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  3.697211 4.574490
## sample estimates:
## difference in location 
##               4.299239
```

#### Ternary plots
For the arid samples

```r
Rock_weathering_filt3_GMPR_Arid_rel <- transform_sample_counts(Rock_weathering_filt3_GMPR_Arid, function(x) x / sum(x) ) # rel abundance
Rock_weathering_filt3_GMPR_Arid_merged <- merge_samples(Rock_weathering_filt3_GMPR_Arid_rel, "Source", fun = mean) # merge by source
Rock_weathering_filt3_GMPR_Arid_merged_rel <- transform_sample_counts(Rock_weathering_filt3_GMPR_Arid_merged, function(x) x / sum(x) ) # rel abundance per source
meandf <- as(otu_table(Rock_weathering_filt3_GMPR_Arid_merged_rel), "matrix")
if (!taxa_are_rows(Rock_weathering_filt3_GMPR_Arid_merged_rel)) { meandf <- t(meandf) }
abundance <- rowSums(meandf) / sum(meandf) * 100

Arid4Ternary <- data.frame(
  meandf,
  Abundance = abundance,
  Phylum = tax_table(Rock_weathering_filt3_GMPR_Arid_merged_rel)[, "Phylum"]
)
# Arid4Ternary <- dplyr::rename(Arid4Ternary, Loess_soil = Loess.soil)

Arid4Ternary$Phylum <-
  factor(Arid4Ternary$Phylum, levels = c(levels(Arid4Ternary$Phylum), 'Rare'))
Arid4Ternary$Phylum[Arid4Ternary$Phylum %in% Rare_phyla]  <- "Rare"
Arid4Ternary$Phylum %<>% 
  factor(., levels = Taxa_rank$Phylum) %>% 
  fct_relevel(., "Rare", after = Inf)

p_ternary_arid <-
  ggtern(data = Arid4Ternary,
         aes(
           x = Loess.soil,
           y = Dust,
           z = Limestone,
           size = Abundance,
           colour = Phylum
         )) +
  geom_point(alpha = 1 / 2) +
  scale_size(
    range = c(1, 5),
    name = "Abundance (%)"
  ) +
  theme_arrownormal() +
    scale_color_manual(values = pal("d3js")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Loess soil") + 
  theme(axis.title = element_blank())
print(p_ternary_arid)
```

![](Rock_weathering_figures/arid ternary-1.svg)<!-- -->

For the hyperarid samples

```r
Rock_weathering_filt3_GMPR_Hyperarid_rel <- transform_sample_counts(Rock_weathering_filt3_GMPR_Hyperarid, function(x) x / sum(x) ) # rel abundance
Rock_weathering_filt3_GMPR_Hyperarid_merged <- merge_samples(Rock_weathering_filt3_GMPR_Hyperarid_rel, "Source", fun = mean) # merge by source
Rock_weathering_filt3_GMPR_Hyperarid_merged_rel <- transform_sample_counts(Rock_weathering_filt3_GMPR_Hyperarid_merged, function(x) x / sum(x) ) # rel abundance per source
meandf <- as(otu_table(Rock_weathering_filt3_GMPR_Hyperarid_merged_rel), "matrix")
if (!taxa_are_rows(Rock_weathering_filt3_GMPR_Hyperarid_merged_rel)) { meandf <- t(meandf) }
abundance <- rowSums(meandf) / sum(meandf) * 100

Hyperarid4Ternary <- data.frame(
  meandf,
  Abundance = abundance,
  Phylum = tax_table(Rock_weathering_filt3_GMPR_Hyperarid_merged_rel)[, "Phylum"]
)

Hyperarid4Ternary$Phylum <-
  factor(Hyperarid4Ternary$Phylum, levels = c(levels(Hyperarid4Ternary$Phylum), 'Rare'))
Hyperarid4Ternary$Phylum[Hyperarid4Ternary$Phylum %in% Rare_phyla]  <- "Rare"
Hyperarid4Ternary$Phylum %<>% 
  factor(., levels = Taxa_rank$Phylum) %>% 
  fct_relevel(., "Rare", after = Inf)

p_ternary_hyperarid <-
  ggtern(data = Hyperarid4Ternary,
         aes(
           x = Loess.soil,
           y = Dust,
           z = Dolomite,
           size = Abundance,
           colour = Phylum
         )) +
  geom_point(alpha = 1 / 2) +
  scale_size(
    range = c(1, 5),
    name = "Abundance (%)"
  ) +
  theme_arrownormal() +
  scale_color_manual(values = pal("d3js")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Loess soil") + 
  theme(axis.title = element_blank())
print(p_ternary_hyperarid)
```

![](Rock_weathering_figures/hyperarid ternary-1.svg)<!-- -->

#### Combine plots
Combine all sequence analysis plots to make Fig. 3

```r
ternary_legend <-
  get_legend(p_ternary_arid)# + theme(legend.direction = "horizontal"))
ord_legend <- get_legend(p_ord)
top_row <-
  plot_grid(
    p_alpha + theme(
      legend.position = "none",
      panel.spacing = unit(0.5, "lines")
    ),
    p_ord + theme(axis.title.y = element_text(vjust = -3)) ,
    labels = c('A', 'B'),
    label_size = 12,
    align = 'v',
    axis = "tl",
    nrow = 1,
    ncol = 2
  )
bottom_l <-
  plot_grid(
    p_taxa_box + theme(legend.position = "none"),
    labels = c('C'),
    label_size = 12,
    ncol = 1
  )
bottom_r <-
  plot_grid(
    p_ternary_arid + 
      theme(legend.position = "none", 
            plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "cm"),
            axis.title = element_blank()),
    p_ternary_hyperarid + 
      theme(legend.position = "none", 
            plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "cm"),
            axis.title = element_blank()),
    labels = c('D'),
    label_size = 12,
    align = 'hv',
    axis = "t",
    # rel_widths = c(1, 1, 0.1),
    scale = c(1.2, 1.2),
    nrow = 2,
    ncol = 1
  )

bottom_rows <- plot_grid(bottom_l, 
                         bottom_r,
                         ternary_legend,
                         align = 'h',
                         axis = "l",
                         scale = c(1, 1, 0.08),
                         rel_widths = c(0.5, 0.35, 0.15),
                         nrow = 1, 
                         ncol = 3)
p_all <- plot_grid(top_row, bottom_rows, align = 'v', axis = 'l', nrow = 2, rel_heights = c(0.43, 0.6)) # aligning vertically along the left axis
print(p_all)
```

![](Rock_weathering_figures/combined plots-1.svg)<!-- -->

### Differential abundance models
Detect differentially abundant OTUs using ALDEx2 [@fernandes_anova-like_2013]

```r
# Rock_weathering_filt3_s <- prune_taxa(names(sort(taxa_sums(Rock_weathering_filt3), TRUE))[1:100], Rock_weathering_filt3)

# run full model 
data2test <- t(otu_table(Rock_weathering_filt3))
comparison <- as.character(unlist(sample_data(Rock_weathering_filt3)[, "Climate.Source"]))
ALDEx_full <- aldex.clr(data2test, comparison, mc.samples = 128, denom = "iqlr", verbose = TRUE, useMC = TRUE) # iqlr for slight assymetry in composition
```

```
## [1] "multicore environment is is OK -- using the BiocParallel package"
## [1] "removed rows with sums equal to zero"
## [1] "computing iqlr centering"
## [1] "data format is OK"
## [1] "dirichlet samples complete"
## [1] "clr transformation complete"
```

```r
ALDEx_full_glm <- aldex.glm(ALDEx_full, comparison, useMC = TRUE) # for more than two conditions
```

```
## [1] "multicore environment is OK -- using the BiocParallel package"
## [1] "running tests for each MC instance:"
## |------------(25%)----------(50%)----------(75%)----------|
```

```r
sig_taxa <- rownames(ALDEx_full_glm)[ALDEx_full_glm$glm.eBH < 0.05] # save names of taxa that are significant under the full model

# Pairwise comparisons
# 
# dolomite - limestone
Rock_weathering_filt3_Rocks <- subset_samples(Rock_weathering_filt3, Uni.Source == "Rock")
ALDEx2plot_Rocks <- CalcALDEx(Rock_weathering_filt3_Rocks, sig_level = 0.1, LFC = 0)
```

```
## [1] "multicore environment is is OK -- using the BiocParallel package"
## [1] "removed rows with sums equal to zero"
## [1] "computing iqlr centering"
## [1] "data format is OK"
## [1] "dirichlet samples complete"
## [1] "clr transformation complete"
## [1] "running tests for each MC instance:"
## |------------(25%)----------(50%)----------(75%)----------|
## [1] "multicore environment is OK -- using the BiocParallel package"
## [1] "sanity check complete"
## [1] "rab.all  complete"
## [1] "rab.win  complete"
## [1] "rab of samples complete"
## [1] "within sample difference calculated"
## [1] "between group difference calculated"
## [1] "group summaries calculated"
## [1] "effect size calculated"
## [1] "summarizing output"
```

```r
GGPlotALDExTax(ALDEx2plot_Rocks) + 
  ggtitle("Hyperarid dolomite vs. Arid limestone")
```

![](Rock_weathering_figures/ALDEx2-1.svg)<!-- -->

```r
# dolomite - soil
Rock_weathering_filt3_DolSoil <- subset_samples(Rock_weathering_filt3, Climate == "Hyperarid" & Source != "Dust")
ALDEx2plot_DolSoil <- CalcALDEx(Rock_weathering_filt3_DolSoil, sig_level = 0.1, LFC = 0)
```

```
## [1] "multicore environment is is OK -- using the BiocParallel package"
## [1] "removed rows with sums equal to zero"
## [1] "computing iqlr centering"
## [1] "data format is OK"
## [1] "dirichlet samples complete"
## [1] "clr transformation complete"
## [1] "running tests for each MC instance:"
## |------------(25%)----------(50%)----------(75%)----------|
## [1] "multicore environment is OK -- using the BiocParallel package"
## [1] "sanity check complete"
## [1] "rab.all  complete"
## [1] "rab.win  complete"
## [1] "rab of samples complete"
## [1] "within sample difference calculated"
## [1] "between group difference calculated"
## [1] "group summaries calculated"
## [1] "effect size calculated"
## [1] "summarizing output"
```

```r
GGPlotALDExTax(ALDEx2plot_DolSoil) + 
  ggtitle("Hyperarid dolomite vs. Hyperarid soil")
```

![](Rock_weathering_figures/ALDEx2-2.svg)<!-- -->

```r
# dolomite - dust
Rock_weathering_filt3_DolDust <- subset_samples(Rock_weathering_filt3, Climate == "Hyperarid" & Source != "Loess soil")
ALDEx2plot_DolDust <- CalcALDEx(Rock_weathering_filt3_DolDust, sig_level = 0.3, LFC = 0)
```

```
## [1] "multicore environment is is OK -- using the BiocParallel package"
## [1] "removed rows with sums equal to zero"
## [1] "computing iqlr centering"
## [1] "data format is OK"
## [1] "dirichlet samples complete"
## [1] "clr transformation complete"
## [1] "running tests for each MC instance:"
## |------------(25%)----------(50%)----------(75%)----------|
## [1] "multicore environment is OK -- using the BiocParallel package"
## [1] "sanity check complete"
## [1] "rab.all  complete"
## [1] "rab.win  complete"
## [1] "rab of samples complete"
## [1] "within sample difference calculated"
## [1] "between group difference calculated"
## [1] "group summaries calculated"
## [1] "effect size calculated"
## [1] "summarizing output"
```

```r
GGPlotALDExTax(ALDEx2plot_DolDust) + 
  ggtitle("Hyperarid dolomite vs. Hyperarid dust")
```

![](Rock_weathering_figures/ALDEx2-3.svg)<!-- -->

```r
# limestone - soil
Rock_weathering_filt3_LimeSoil <- subset_samples(Rock_weathering_filt3, Climate == "Arid" & Source != "Dust")
ALDEx2plot_LimeSoil <- CalcALDEx(Rock_weathering_filt3_LimeSoil, sig_level = 0.1, LFC = 0)
```

```
## [1] "multicore environment is is OK -- using the BiocParallel package"
## [1] "removed rows with sums equal to zero"
## [1] "computing iqlr centering"
## [1] "data format is OK"
## [1] "dirichlet samples complete"
## [1] "clr transformation complete"
## [1] "running tests for each MC instance:"
## |------------(25%)----------(50%)----------(75%)----------|
## [1] "multicore environment is OK -- using the BiocParallel package"
## [1] "sanity check complete"
## [1] "rab.all  complete"
## [1] "rab.win  complete"
## [1] "rab of samples complete"
## [1] "within sample difference calculated"
## [1] "between group difference calculated"
## [1] "group summaries calculated"
## [1] "effect size calculated"
## [1] "summarizing output"
```

```r
GGPlotALDExTax(ALDEx2plot_LimeSoil) + 
  ggtitle("Arid limestone vs. Arid soil")
```

![](Rock_weathering_figures/ALDEx2-4.svg)<!-- -->

```r
# limestone - dust
Rock_weathering_filt3_LimeDust <- subset_samples(Rock_weathering_filt3, Climate == "Arid" & Source != "Loess soil")
ALDEx2plot_LimeDust <- CalcALDEx(Rock_weathering_filt3_LimeDust, sig_level = 0.3, LFC = 0)
```

```
## [1] "multicore environment is is OK -- using the BiocParallel package"
## [1] "removed rows with sums equal to zero"
## [1] "computing iqlr centering"
## [1] "data format is OK"
## [1] "dirichlet samples complete"
## [1] "clr transformation complete"
## [1] "running tests for each MC instance:"
## |------------(25%)----------(50%)----------(75%)----------|
## [1] "multicore environment is OK -- using the BiocParallel package"
## [1] "sanity check complete"
## [1] "rab.all  complete"
## [1] "rab.win  complete"
## [1] "rab of samples complete"
## [1] "within sample difference calculated"
## [1] "between group difference calculated"
## [1] "group summaries calculated"
## [1] "effect size calculated"
## [1] "summarizing output"
```

```r
GGPlotALDExTax(ALDEx2plot_LimeDust) + 
  ggtitle("Arid limestone vs. Arid dust")
```

![](Rock_weathering_figures/ALDEx2-5.svg)<!-- -->

```r
ALDEx2plot_Rocks %<>% cbind(., Var1 = "Dolomite", Var2 = "Limestone")
ALDEx2plot_DolSoil %<>% cbind(., Var1 = "Dolomite", Var2 = "Loess soil")
ALDEx2plot_DolDust %<>% cbind(., Var1 = "Dolomite", Var2 = "Dust")
ALDEx2plot_LimeSoil %<>% cbind(., Var1 = "Limestone", Var2 = "Loess soil")
ALDEx2plot_LimeDust %<>% cbind(., Var1 = "Limestone", Var2 = "Dust")

ALDEx2plot_all <- bind_rows(ALDEx2plot_Rocks, ALDEx2plot_DolSoil, ALDEx2plot_DolDust, ALDEx2plot_LimeSoil, ALDEx2plot_LimeDust)
ALDEx2plot_all$Var2 %<>%
    factor() %>%  # Taxa_rank is calcuted for the taxa box plots
    fct_relevel(., "Limestone")

# paste0(percent(sum(ALDEx2plot_Rocks$effect > 0 & ALDEx2plot_Rocks$Significance == "Pass")/nrow(ALDEx2plot_Rocks)), "/", percent(sum(ALDEx2plot_Rocks$effect < 0 & ALDEx2plot_Rocks$Significance == "Pass")/nrow(ALDEx2plot_Rocks)))

Labels <- c(
  paste0("⬆", sum(ALDEx2plot_Rocks$effect > 0 & ALDEx2plot_Rocks$Significance == "Pass"), " ⬇", sum(ALDEx2plot_Rocks$effect < 0 & ALDEx2plot_Rocks$Significance == "Pass"), " (", nrow(ALDEx2plot_Rocks), ")"),
  paste0("⬆", sum(ALDEx2plot_DolSoil$effect > 0 & ALDEx2plot_DolSoil$Significance == "Pass"), " ⬇", sum(ALDEx2plot_DolSoil$effect < 0 & ALDEx2plot_DolSoil$Significance == "Pass"), " (", nrow(ALDEx2plot_DolSoil), ")"),
  paste0("⬆", sum(ALDEx2plot_DolDust$effect > 0 & ALDEx2plot_DolDust$Significance == "Pass"), " ⬇", sum(ALDEx2plot_DolDust$effect < 0 & ALDEx2plot_DolDust$Significance == "Pass"), " (", nrow(ALDEx2plot_DolDust), ")"),
  paste0("⬆", sum(ALDEx2plot_LimeSoil$effect > 0 & ALDEx2plot_LimeSoil$Significance == "Pass"), " ⬇", sum(ALDEx2plot_LimeSoil$effect < 0 & ALDEx2plot_LimeSoil$Significance == "Pass"), " (", nrow(ALDEx2plot_LimeSoil), ")"),
  paste0("⬆", sum(ALDEx2plot_LimeDust$effect > 0 & ALDEx2plot_LimeDust$Significance == "Pass"), " ⬇", sum(ALDEx2plot_LimeDust$effect < 0 & ALDEx2plot_LimeDust$Significance == "Pass"), " (", nrow(ALDEx2plot_LimeDust), ")")
)
Label_text <- bind_cols(
  unique(ALDEx2plot_all[c("Var1", "Var2")]),
  Label = Labels
  )
```

```r
p_aldex2_all <- GGPlotALDExTax(ALDEx2plot_all) +
  facet_grid(Var2 ~ Var1, scales = "free_y") +
  # theme(strip.background = element_blank(), strip.placement = "outside") +
  geom_text(
    data    = Label_text,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  ) 
print(p_aldex2_all)
```

![](Rock_weathering_figures/ALDEx2_combined_plot-1.svg)<!-- -->

## Other plots
Other plots in the paper which are not based on sequence data
### Isotopes profile

```r
Isotopes <-
  read_csv(
    "Data/Isotopes_data.csv"
  )

Isotopes %<>% 
  mutate(Mean.Arid = (`Limestone Shivta Fm. NWSH1` + `Limestone Shivta Fm. NWSH2`) / 2)
Isotopes %<>% 
  mutate(Mean.Hyperarid = (`Dolomite Gerofit Fm.UVSL5` + `Dolomite Gerofit Fm.UVSL6` ) / 2)

Isotopes2plot <- data.frame(
  Rock = factor(c(rep("Limestone", 10), rep("Dolomite", 10)), 
                levels = c("Limestone", "Dolomite")),
  Depth = rep(Isotopes$`Depth (mm)`, 2),
  Isotope = rep(Isotopes$Isotope, 2),
  min = c(
    pmin(
      Isotopes$`Limestone Shivta Fm. NWSH1`,
      Isotopes$`Limestone Shivta Fm. NWSH2`
    ),
    pmin(
      Isotopes$`Dolomite Gerofit Fm.UVSL5`,
      Isotopes$`Dolomite Gerofit Fm.UVSL6`
    )
  ),
  max = c(
    pmax(
      Isotopes$`Limestone Shivta Fm. NWSH1`,
      Isotopes$`Limestone Shivta Fm. NWSH2`
    ),
    pmax(
      Isotopes$`Dolomite Gerofit Fm.UVSL5`,
      Isotopes$`Dolomite Gerofit Fm.UVSL6`
    )
  ),
  mean = c(Isotopes$Mean.Arid, Isotopes$Mean.Hyperarid)
)

p_isotopes <-
  ggplot(Isotopes2plot, aes(y = mean, x = Depth, colour = Isotope)) +
  geom_point(size = 4, alpha = 1 / 2) +
  geom_errorbar(aes(ymin = min, ymax = max), alpha = 1/2, width = 0.2) +
  geom_line(alpha = 1 / 2) +
  coord_flip() +
  theme_cowplot(font_size = 18, font_family = f_name) +
  background_grid(major = "xy",
                  minor = "none") +
  scale_x_reverse(limits = c(4.1, -0.1), expand = c(0.01, 0.01)) +
  # scale_x_continuous(limits = c(0, 50), expand = c(0.01, 0.01)) +
  facet_grid(Rock ~ . , scales = "free_x", labeller = label_parsed) +
  scale_color_manual(values =  pom4[c(2,1)],
                     labels = c(expression(paste(delta ^ {13}, "C")),
                                expression(paste(delta ^ {18}, "O")))) +
  ylab(expression(paste(delta ^ {13}, "C / ",
                        delta ^ {18}, "O", " (", "\u2030", "VPDB",")"
  )))

p_isotopes <- plot_grid(p_isotopes, labels = "b", label_size = 20)
print(p_isotopes)
```

![](Rock_weathering_figures/isotopes-1.svg)<!-- -->

### Drying experiment

```r
read_csv("Data/Drying_data_full.csv") ->
  Drying_long
Drying_long$Rock %<>% fct_relevel(., "Limestone")
Drying_long$BRC %<>% 
  fct_relevel(., "Present")
Drying_long$Sample <- with(Drying_long, paste(Rock, BRC))

Drying_mods <- tibble(Sample = character(), Intercept = numeric(), b = numeric(), a = numeric(), P = numeric(), R2 = numeric())
mods <- list()
j <- 1
for (i in unique(Drying_long$Sample)) {
  data2model <- Drying_long[Drying_long$Sample == i, ]
  colnames(data2model) <- c("Time", "Replicate", "BRC", "Rock", "RWC", "Sample")
  (mod <- lme(RWC ~ poly(Time, 2, raw = TRUE), random = ~0 + Time|Replicate, data = data2model))
  intervals(mod)
  # mod <- lm(`Residual water content (%)` ~ sqrt(1/(`Time (h)` + 1)), data = data2model)
  mods[[j]] <- mod
  Drying_mods[j, "Sample"] <- i
  Drying_mods[j, "Intercept"] <- mod$coefficients$fixed[1]
  Drying_mods[j, "b"] <- mod$coefficients$fixed[2]
  Drying_mods[j, "a"] <- mod$coefficients$fixed[3]
  Drying_mods[j, "P"] <- anova(mod)$`p-value`[2]
  Drying_mods[j, "R2"] <- r.squaredGLMM(mod)[, "R2c"]
  j <- j + 1
}

Drying_mods %>% 
  kable(., digits = c(1, 1, 2, 2, 3, 2)) %>%
  kable_styling(bootstrap_options = c("hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:right;"> Intercept </th>
   <th style="text-align:right;"> b </th>
   <th style="text-align:right;"> a </th>
   <th style="text-align:right;"> P </th>
   <th style="text-align:right;"> R2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Dolomite Present </td>
   <td style="text-align:right;"> 97.2 </td>
   <td style="text-align:right;"> -0.85 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.002 </td>
   <td style="text-align:right;"> 0.86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Dolomite Removed </td>
   <td style="text-align:right;"> 91.9 </td>
   <td style="text-align:right;"> -2.42 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.004 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Limestone Present </td>
   <td style="text-align:right;"> 97.4 </td>
   <td style="text-align:right;"> -0.99 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Limestone Removed </td>
   <td style="text-align:right;"> 92.5 </td>
   <td style="text-align:right;"> -3.83 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.92 </td>
  </tr>
</tbody>
</table>

```r
# comapre with and without crust
data2model <- Drying_long
colnames(data2model) <- c("Time", "Replicate", "BRC", "Rock", "RWC", "Sample")
mod_all <- lme(RWC ~ poly(Time, 2, raw = TRUE), random = ~0 + Time|Replicate, data = data2model)
mod_treatment <- lme(RWC ~ poly(Time, 2, raw = TRUE) * BRC, random = ~0 + Time|Replicate, data = data2model)
anova(mod_all, mod_treatment)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> call </th>
   <th style="text-align:right;"> Model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AIC </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> logLik </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;"> L.Ratio </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mod_all </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE), data = data2model,     random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1108.742 </td>
   <td style="text-align:right;"> 1122.884 </td>
   <td style="text-align:right;"> -549.3711 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mod_treatment </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE) * BRC, data = data2model,     random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1024.040 </td>
   <td style="text-align:right;"> 1046.472 </td>
   <td style="text-align:right;"> -504.0198 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 90.70266 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

</div>

```r
# comapre limestone vs dolomite - with crust
mod_all <- lme(RWC ~ poly(Time, 2, raw = TRUE), random = ~0 + Time|Replicate, data = data2model[data2model$BRC == "Present", ])
mod_treatment <- lme(RWC ~ poly(Time, 2, raw = TRUE) * Rock, random = ~0 + Time|Replicate, data = data2model[data2model$BRC == "Present", ])
anova(mod_all, mod_treatment)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> call </th>
   <th style="text-align:right;"> Model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AIC </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> logLik </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;"> L.Ratio </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mod_all </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE), data = data2model[data2model$BRC ==     &quot;Present&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 409.0114 </td>
   <td style="text-align:right;"> 419.5658 </td>
   <td style="text-align:right;"> -199.5057 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mod_treatment </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE) * Rock, data = data2model[data2model$BRC ==     &quot;Present&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 409.0676 </td>
   <td style="text-align:right;"> 425.5512 </td>
   <td style="text-align:right;"> -196.5338 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 5.943792 </td>
   <td style="text-align:right;"> 0.1143771 </td>
  </tr>
</tbody>
</table>

</div>

```r
# comapre limestone vs dolomite - without crust
mod_all <- lme(RWC ~ poly(Time, 2, raw = TRUE), random = ~0 + Time|Replicate, data = data2model[data2model$BRC == "Removed", ])
mod_treatment <- lme(RWC ~ poly(Time, 2, raw = TRUE) * Rock, random = ~0 + Time|Replicate, data = data2model[data2model$BRC == "Removed", ])
anova(mod_all, mod_treatment)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> call </th>
   <th style="text-align:right;"> Model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AIC </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> logLik </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;"> L.Ratio </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mod_all </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE), data = data2model[data2model$BRC ==     &quot;Removed&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 530.3667 </td>
   <td style="text-align:right;"> 540.9211 </td>
   <td style="text-align:right;"> -260.1834 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mod_treatment </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE) * Rock, data = data2model[data2model$BRC ==     &quot;Removed&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 525.9703 </td>
   <td style="text-align:right;"> 542.4538 </td>
   <td style="text-align:right;"> -254.9851 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 10.39642 </td>
   <td style="text-align:right;"> 0.0154802 </td>
  </tr>
</tbody>
</table>

</div>

```r
# comapre with and without crust - limestone
mod_all <- lme(RWC ~ poly(Time, 2, raw = TRUE), random = ~0 + Time|Replicate, data = data2model[data2model$Rock == "Limestone", ])
mod_treatment <- lme(RWC ~ poly(Time, 2, raw = TRUE) * BRC, random = ~0 + Time|Replicate, data = data2model[data2model$Rock == "Limestone", ])
anova(mod_all, mod_treatment)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> call </th>
   <th style="text-align:right;"> Model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AIC </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> logLik </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;"> L.Ratio </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mod_all </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE), data = data2model[data2model$Rock ==     &quot;Limestone&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 562.9452 </td>
   <td style="text-align:right;"> 573.4996 </td>
   <td style="text-align:right;"> -276.4726 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mod_treatment </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE) * BRC, data = data2model[data2model$Rock ==     &quot;Limestone&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 497.8492 </td>
   <td style="text-align:right;"> 514.3328 </td>
   <td style="text-align:right;"> -240.9246 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 71.09599 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

</div>

```r
mod_all <- lme(RWC ~ poly(Time, 2, raw = TRUE), random = ~0 + Time|Replicate, data = data2model[data2model$Rock == "Dolomite", ])
mod_treatment <- lme(RWC ~ poly(Time, 2, raw = TRUE) * BRC, random = ~0 + Time|Replicate, data = data2model[data2model$Rock == "Dolomite", ])
anova(mod_all, mod_treatment)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> call </th>
   <th style="text-align:right;"> Model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AIC </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> logLik </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;"> L.Ratio </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mod_all </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE), data = data2model[data2model$Rock ==     &quot;Dolomite&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 555.4483 </td>
   <td style="text-align:right;"> 566.0027 </td>
   <td style="text-align:right;"> -272.7241 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mod_treatment </td>
   <td style="text-align:left;"> lme.formula(fixed = RWC ~ poly(Time, 2, raw = TRUE) * BRC, data = data2model[data2model$Rock ==     &quot;Dolomite&quot;, ], random = ~0 + Time | Replicate) </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 529.4140 </td>
   <td style="text-align:right;"> 545.8975 </td>
   <td style="text-align:right;"> -256.7070 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 32.03428 </td>
   <td style="text-align:right;"> 5e-07 </td>
  </tr>
</tbody>
</table>

</div>

```r
p_drying <-
  ggplot(
    Drying_long,
    aes(
      x = `Time (h)`,
      y = `Residual water content (%)`,
      colour = Rock,
      # fill = Rock,
      shape = BRC,
      linetype = BRC
    )
  ) +
  geom_point(size = 2, alpha = 2/3) +
  # geom_smooth(method = "lm", se = FALSE, alpha = 1/2, formula = (y ~ sqrt(1/(x+1)))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 1/3, formula = (y ~ poly(x, 2)), size = 1) +
  # geom_line(alpha = 1/2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0.01, 0.01)) +
  # scale_fill_manual(values = pom4) +
  scale_color_manual(values = pom4)
print(p_drying)
```

![](Rock_weathering_figures/drying-1.svg)<!-- -->


```r
devtools::session_info()
```

```
## ─ Session info ─────────────────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 3.4.4 (2018-03-15)
##  os       KDE neon User Edition 5.14  
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language en_GB                       
##  collate  en_DK.UTF-8                 
##  ctype    en_DK.UTF-8                 
##  tz       Europe/Vienna               
##  date     2019-01-04                  
## 
## ─ Packages ─────────────────────────────────────────────────────────────────────────────
##  package              * version   date       lib source                                  
##  abind                  1.4-5     2016-07-21 [1] CRAN (R 3.4.4)                          
##  acepack                1.4.1     2016-10-29 [1] CRAN (R 3.5.1)                          
##  ade4                   1.7-13    2018-08-31 [1] CRAN (R 3.5.1)                          
##  affy                   1.58.0    2018-09-25 [1] Bioconductor                            
##  affyio                 1.50.0    2018-09-25 [1] Bioconductor                            
##  agricolae            * 1.2-8     2017-09-12 [1] CRAN (R 3.4.4)                          
##  ALDEx2               * 1.10.0    2019-01-04 [1] Bioconductor                            
##  AlgDesign              1.1-7.3   2014-10-15 [1] CRAN (R 3.5.1)                          
##  ape                    5.2       2018-09-24 [1] CRAN (R 3.4.4)                          
##  artyfarty            * 0.0.1     2018-09-25 [1] Github (datarootsio/artyfarty@e2b3804)  
##  assertthat             0.2.0     2017-04-11 [1] CRAN (R 3.5.1)                          
##  backports              1.1.3     2018-12-14 [1] CRAN (R 3.4.4)                          
##  base64enc              0.1-3     2015-07-28 [1] CRAN (R 3.5.1)                          
##  bayesm                 3.1-1     2018-12-21 [1] CRAN (R 3.4.4)                          
##  BiasedUrn              1.07      2015-12-28 [1] CRAN (R 3.4.4)                          
##  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.1)                          
##  bindrcpp             * 0.2.2     2018-03-29 [1] CRAN (R 3.4.4)                          
##  Biobase              * 2.38.0    2019-01-04 [1] Bioconductor                            
##  BiocGenerics         * 0.24.0    2019-01-04 [1] Bioconductor                            
##  BiocInstaller          1.28.0    2019-01-04 [1] Bioconductor                            
##  BiocParallel           1.12.0    2019-01-04 [1] Bioconductor                            
##  BiodiversityR        * 2.10-1    2018-07-14 [1] CRAN (R 3.4.4)                          
##  biomformat             1.8.0     2018-09-26 [1] Bioconductor                            
##  Biostrings             2.48.0    2018-09-25 [1] Bioconductor                            
##  bitops                 1.0-6     2013-08-17 [1] CRAN (R 3.5.1)                          
##  boot                   1.3-20    2017-07-30 [4] CRAN (R 3.4.2)                          
##  brew                   1.0-6     2011-04-13 [1] CRAN (R 3.5.1)                          
##  broom                  0.5.1     2018-12-05 [1] CRAN (R 3.4.4)                          
##  callr                  3.1.1     2018-12-21 [1] CRAN (R 3.4.4)                          
##  car                  * 3.0-2     2018-08-23 [1] CRAN (R 3.4.4)                          
##  carData              * 3.0-2     2018-09-30 [1] CRAN (R 3.4.4)                          
##  cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.5.1)                          
##  checkmate              1.8.5     2017-10-24 [1] CRAN (R 3.4.4)                          
##  class                  7.3-15    2019-01-01 [1] CRAN (R 3.4.4)                          
##  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)                          
##  cluster                2.0.7-1   2018-04-09 [1] CRAN (R 3.5.1)                          
##  coda                   0.19-2    2018-10-08 [1] CRAN (R 3.4.4)                          
##  codetools              0.2-16    2018-12-24 [4] CRAN (R 3.4.4)                          
##  coin                   1.2-2     2017-11-28 [1] CRAN (R 3.5.1)                          
##  colorspace             1.3-2     2016-12-14 [1] CRAN (R 3.5.1)                          
##  combinat               0.0-8     2012-10-29 [1] CRAN (R 3.5.1)                          
##  compositions           1.40-2    2018-06-14 [1] CRAN (R 3.5.1)                          
##  cowplot              * 0.9.3     2018-07-15 [1] CRAN (R 3.5.1)                          
##  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.1)                          
##  curl                   3.2       2018-03-28 [1] CRAN (R 3.5.1)                          
##  data.table             1.11.8    2018-09-30 [1] CRAN (R 3.4.4)                          
##  data.tree              0.7.8     2018-09-24 [1] CRAN (R 3.4.4)                          
##  DelayedArray           0.6.6     2018-09-11 [1] Bioconductor                            
##  deldir                 0.1-15    2018-04-01 [1] CRAN (R 3.5.1)                          
##  DEoptimR               1.0-8     2016-11-19 [1] CRAN (R 3.5.1)                          
##  desc                   1.2.0     2018-05-01 [1] CRAN (R 3.5.1)                          
##  DescTools              0.99.26   2018-11-13 [1] CRAN (R 3.4.4)                          
##  devtools             * 2.0.1     2018-10-26 [1] CRAN (R 3.4.4)                          
##  DiagrammeR             1.0.0     2018-03-01 [1] CRAN (R 3.4.4)                          
##  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.4.4)                          
##  diptest                0.75-7    2016-12-05 [1] CRAN (R 3.4.4)                          
##  doParallel           * 1.0.14    2018-09-24 [1] CRAN (R 3.5.1)                          
##  downloader             0.4       2015-07-09 [1] CRAN (R 3.5.1)                          
##  dplyr                * 0.7.8     2018-11-10 [1] CRAN (R 3.4.4)                          
##  e1071                  1.7-0     2018-07-28 [1] CRAN (R 3.5.1)                          
##  effects                4.0-3     2018-08-19 [1] CRAN (R 3.5.1)                          
##  EMT                    1.1       2013-01-29 [1] CRAN (R 3.5.1)                          
##  energy                 1.7-5     2018-08-11 [1] CRAN (R 3.4.4)                          
##  evaluate               0.12      2018-10-09 [1] CRAN (R 3.4.4)                          
##  expm                   0.999-3   2018-09-22 [1] CRAN (R 3.5.1)                          
##  extrafont            * 0.17      2014-12-08 [1] CRAN (R 3.5.1)                          
##  extrafontdb            1.0       2012-06-11 [1] CRAN (R 3.5.1)                          
##  forcats              * 0.3.0     2018-02-19 [1] CRAN (R 3.4.4)                          
##  foreach              * 1.4.4     2017-12-12 [1] CRAN (R 3.5.1)                          
##  foreign                0.8-71    2018-07-20 [4] CRAN (R 3.4.4)                          
##  Formula                1.2-3     2018-05-03 [1] CRAN (R 3.5.1)                          
##  fs                     1.2.6     2018-08-23 [1] CRAN (R 3.4.4)                          
##  gdata                  2.18.0    2017-06-06 [1] CRAN (R 3.5.1)                          
##  gdtools              * 0.1.7     2018-02-27 [1] CRAN (R 3.4.4)                          
##  generics               0.0.2     2018-11-29 [1] CRAN (R 3.4.4)                          
##  GenomeInfoDb           1.16.0    2018-09-25 [1] Bioconductor                            
##  GenomeInfoDbData       1.1.0     2018-09-25 [1] Bioconductor                            
##  GenomicRanges          1.32.7    2018-09-20 [1] Bioconductor                            
##  GGally                 1.4.0     2018-05-17 [1] CRAN (R 3.4.4)                          
##  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.4.4)                          
##  ggpomological        * 0.1.2     2019-01-04 [1] Github (gadenbuie/ggpomological@5d0c335)
##  ggrepel                0.8.0     2018-05-09 [1] CRAN (R 3.4.4)                          
##  ggridges               0.5.1     2018-09-27 [1] CRAN (R 3.4.4)                          
##  ggtern               * 3.1.0     2018-12-19 [1] CRAN (R 3.4.4)                          
##  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)                          
##  gmodels                2.18.1    2018-06-25 [1] CRAN (R 3.5.1)                          
##  GPArotation            2014.11-1 2014-11-25 [1] CRAN (R 3.4.4)                          
##  gridExtra              2.3       2017-09-09 [1] CRAN (R 3.4.4)                          
##  gtable                 0.2.0     2016-02-26 [1] CRAN (R 3.4.4)                          
##  gtools                 3.8.1     2018-06-26 [1] CRAN (R 3.5.1)                          
##  haven                  2.0.0     2018-11-22 [1] CRAN (R 3.4.4)                          
##  hexbin               * 1.27.2    2018-01-15 [1] CRAN (R 3.4.4)                          
##  highr                  0.7       2018-06-09 [1] CRAN (R 3.5.1)                          
##  Hmisc                  4.1-1     2018-01-03 [1] CRAN (R 3.5.1)                          
##  hms                    0.4.2     2018-03-10 [1] CRAN (R 3.5.1)                          
##  htmlTable              1.13      2019-01-02 [1] CRAN (R 3.4.4)                          
##  htmltools              0.3.6     2017-04-28 [1] CRAN (R 3.4.4)                          
##  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.4.4)                          
##  httpuv                 1.4.5.1   2018-12-18 [1] CRAN (R 3.4.4)                          
##  httr                   1.4.0     2018-12-11 [1] CRAN (R 3.4.4)                          
##  igraph                 1.2.2     2018-07-27 [1] CRAN (R 3.5.1)                          
##  influenceR             0.1.0     2015-09-03 [1] CRAN (R 3.5.1)                          
##  IRanges                2.14.12   2018-09-20 [1] Bioconductor                            
##  iterators            * 1.0.10    2018-07-13 [1] CRAN (R 3.5.1)                          
##  jsonlite               1.6       2018-12-07 [1] CRAN (R 3.4.4)                          
##  kableExtra           * 0.9.0     2018-05-21 [1] CRAN (R 3.5.1)                          
##  klaR                   0.6-14    2018-03-19 [1] CRAN (R 3.5.1)                          
##  knitr                * 1.21      2018-12-10 [1] CRAN (R 3.4.4)                          
##  labeling               0.3       2014-08-23 [1] CRAN (R 3.5.1)                          
##  later                  0.7.5     2018-09-18 [1] CRAN (R 3.4.4)                          
##  latex2exp              0.4.0     2015-11-30 [1] CRAN (R 3.5.1)                          
##  lattice              * 0.20-38   2018-11-04 [1] CRAN (R 3.4.4)                          
##  latticeExtra           0.6-28    2016-02-09 [1] CRAN (R 3.5.1)                          
##  lavaan                 0.6-3     2018-09-22 [1] CRAN (R 3.4.4)                          
##  lazyeval               0.2.1     2017-10-29 [1] CRAN (R 3.5.1)                          
##  LearnBayes             2.15.1    2018-03-18 [1] CRAN (R 3.5.1)                          
##  limma                  3.36.5    2018-09-20 [1] Bioconductor                            
##  lme4                   1.1-19    2018-11-10 [1] CRAN (R 3.4.4)                          
##  lmtest                 0.9-36    2018-04-04 [1] CRAN (R 3.4.4)                          
##  lubridate              1.7.4     2018-04-11 [1] CRAN (R 3.4.4)                          
##  magrittr             * 1.5       2014-11-22 [1] CRAN (R 3.5.1)                          
##  manipulate             1.0.1     2014-12-24 [1] CRAN (R 3.5.1)                          
##  MASS                 * 7.3-51.1  2018-11-01 [1] CRAN (R 3.4.4)                          
##  Matrix                 1.2-15    2018-11-01 [1] CRAN (R 3.4.4)                          
##  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.4.4)                          
##  MBESS                  4.4.3     2018-01-10 [1] CRAN (R 3.4.4)                          
##  memoise                1.1.0     2017-04-21 [1] CRAN (R 3.5.1)                          
##  mgcv                   1.8-26    2018-11-21 [4] CRAN (R 3.4.4)                          
##  mime                   0.6       2018-10-05 [1] CRAN (R 3.4.4)                          
##  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 3.5.1)                          
##  minpack.lm             1.2-1     2016-11-20 [1] CRAN (R 3.4.4)                          
##  minqa                  1.2.4     2014-10-09 [1] CRAN (R 3.4.4)                          
##  mnormt                 1.5-5     2016-10-15 [1] CRAN (R 3.4.4)                          
##  modelr                 0.1.2     2018-05-11 [1] CRAN (R 3.5.1)                          
##  modeltools             0.2-22    2018-07-16 [1] CRAN (R 3.5.1)                          
##  multcomp               1.4-8     2017-11-08 [1] CRAN (R 3.5.1)                          
##  multcompView           0.1-7     2015-07-31 [1] CRAN (R 3.4.4)                          
##  multtest               2.36.0    2018-09-26 [1] Bioconductor                            
##  MuMIn                * 1.42.1    2018-07-23 [1] CRAN (R 3.4.4)                          
##  munsell                0.5.0     2018-06-12 [1] CRAN (R 3.5.1)                          
##  mvtnorm                1.0-8     2018-05-31 [1] CRAN (R 3.5.1)                          
##  nlme                 * 3.1-137   2018-04-07 [1] CRAN (R 3.4.4)                          
##  nloptr                 1.2.1     2018-10-03 [1] CRAN (R 3.4.4)                          
##  nnet                   7.3-12    2016-02-02 [1] CRAN (R 3.5.1)                          
##  nortest                1.0-4     2015-07-30 [1] CRAN (R 3.5.1)                          
##  openxlsx               4.1.0     2018-05-26 [1] CRAN (R 3.4.4)                          
##  pander                 0.6.3     2018-11-06 [1] CRAN (R 3.4.4)                          
##  pbivnorm               0.6.0     2015-01-23 [1] CRAN (R 3.4.4)                          
##  permute              * 0.9-4     2016-09-09 [1] CRAN (R 3.4.4)                          
##  phyloseq             * 1.22.3    2019-01-04 [1] Bioconductor                            
##  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.4.4)                          
##  pkgbuild               1.0.2     2018-10-16 [1] CRAN (R 3.4.4)                          
##  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)                          
##  pkgload                1.0.2     2018-10-29 [1] CRAN (R 3.4.4)                          
##  plyr                   1.8.4     2016-06-08 [1] CRAN (R 3.4.4)                          
##  preprocessCore         1.42.0    2018-09-25 [1] Bioconductor                            
##  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.1)                          
##  processx               3.2.1     2018-12-05 [1] CRAN (R 3.4.4)                          
##  promises               1.0.1     2018-04-13 [1] CRAN (R 3.4.4)                          
##  proto                  1.0.0     2016-10-29 [1] CRAN (R 3.5.1)                          
##  ps                     1.3.0     2018-12-21 [1] CRAN (R 3.4.4)                          
##  psych                  1.8.10    2018-10-31 [1] CRAN (R 3.4.4)                          
##  purrr                * 0.2.5     2018-05-29 [1] CRAN (R 3.4.4)                          
##  pwr                    1.2-2     2018-03-03 [1] CRAN (R 3.4.4)                          
##  questionr              0.7.0     2018-11-26 [1] CRAN (R 3.4.4)                          
##  R6                     2.3.0     2018-10-04 [1] CRAN (R 3.4.4)                          
##  Rcmdr                  2.5-1     2018-09-11 [1] CRAN (R 3.5.1)                          
##  RcmdrMisc              2.5-1     2018-09-10 [1] CRAN (R 3.5.1)                          
##  RColorBrewer           1.1-2     2014-12-07 [1] CRAN (R 3.4.4)                          
##  rcompanion           * 2.0.10    2019-01-03 [1] CRAN (R 3.4.4)                          
##  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.4.4)                          
##  RCurl                  1.95-4.11 2018-07-15 [1] CRAN (R 3.5.1)                          
##  readr                * 1.3.1     2018-12-21 [1] CRAN (R 3.4.4)                          
##  readxl                 1.2.0     2018-12-19 [1] CRAN (R 3.4.4)                          
##  relimp                 1.0-5     2016-03-30 [1] CRAN (R 3.5.1)                          
##  remotes                2.0.2     2018-10-30 [1] CRAN (R 3.4.4)                          
##  reshape                0.8.8     2018-10-23 [1] CRAN (R 3.4.4)                          
##  reshape2               1.4.3     2017-12-11 [1] CRAN (R 3.4.4)                          
##  rgexf                  0.15.3    2015-03-24 [1] CRAN (R 3.5.1)                          
##  rhdf5                  2.24.0    2018-09-25 [1] Bioconductor                            
##  Rhdf5lib               1.2.1     2018-09-25 [1] Bioconductor                            
##  rio                    0.5.16    2018-11-26 [1] CRAN (R 3.4.4)                          
##  rlang                  0.3.0.1   2018-10-25 [1] CRAN (R 3.4.4)                          
##  rmarkdown            * 1.11      2018-12-08 [1] CRAN (R 3.4.4)                          
##  robustbase             0.93-3    2018-09-21 [1] CRAN (R 3.4.4)                          
##  Rook                   1.1-1     2014-10-20 [1] CRAN (R 3.5.1)                          
##  rpart                  4.1-13    2018-02-23 [1] CRAN (R 3.4.4)                          
##  rprojroot              1.3-2     2018-01-03 [1] CRAN (R 3.5.1)                          
##  rstudioapi             0.8       2018-10-02 [1] CRAN (R 3.4.4)                          
##  Rttf2pt1               1.3.7     2018-06-29 [1] CRAN (R 3.5.1)                          
##  rvest                  0.3.2     2016-06-17 [1] CRAN (R 3.5.1)                          
##  S4Vectors              0.18.3    2018-09-25 [1] Bioconductor                            
##  sandwich               2.5-0     2018-08-17 [1] CRAN (R 3.5.1)                          
##  scales               * 1.0.0     2018-08-09 [1] CRAN (R 3.4.4)                          
##  SCRT                   1.2.2     2018-03-07 [1] CRAN (R 3.4.4)                          
##  sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.4.4)                          
##  shiny                  1.2.0     2018-11-02 [1] CRAN (R 3.4.4)                          
##  sp                     1.3-1     2018-06-05 [1] CRAN (R 3.4.4)                          
##  spData                 0.2.9.6   2018-12-03 [1] CRAN (R 3.4.4)                          
##  spdep                  0.8-1     2018-11-21 [1] CRAN (R 3.4.4)                          
##  stringi                1.2.4     2018-07-20 [1] CRAN (R 3.4.4)                          
##  stringr              * 1.3.1     2018-05-10 [1] CRAN (R 3.5.1)                          
##  SummarizedExperiment   1.10.1    2018-09-26 [1] Bioconductor                            
##  SuppDists              1.1-9.4   2016-09-23 [1] CRAN (R 3.4.4)                          
##  survey                 3.35      2018-12-17 [1] CRAN (R 3.4.4)                          
##  survival               2.43-3    2018-11-26 [1] CRAN (R 3.4.4)                          
##  svglite              * 1.2.1     2017-09-11 [1] CRAN (R 3.4.4)                          
##  tcltk2                 1.2-11    2014-12-20 [1] CRAN (R 3.5.1)                          
##  tensorA                0.36.1    2018-07-29 [1] CRAN (R 3.5.1)                          
##  testthat               2.0.1     2018-10-13 [1] CRAN (R 3.4.4)                          
##  TH.data                1.0-9     2018-07-10 [1] CRAN (R 3.5.1)                          
##  tibble               * 1.4.2     2018-01-22 [1] CRAN (R 3.4.4)                          
##  tidyr                * 0.8.2     2018-10-28 [1] CRAN (R 3.4.4)                          
##  tidyselect             0.2.5     2018-10-11 [1] CRAN (R 3.4.4)                          
##  tidyverse            * 1.2.1     2017-11-14 [1] CRAN (R 3.4.4)                          
##  ufs                    0.0.1     2018-08-02 [1] CRAN (R 3.4.4)                          
##  userfriendlyscience  * 0.7.2     2018-09-24 [1] CRAN (R 3.4.4)                          
##  usethis              * 1.4.0     2018-08-14 [1] CRAN (R 3.4.4)                          
##  vegan                * 2.5-3     2018-10-25 [1] CRAN (R 3.4.4)                          
##  viridis                0.5.1     2018-03-29 [1] CRAN (R 3.4.4)                          
##  viridisLite            0.3.0     2018-02-01 [1] CRAN (R 3.5.1)                          
##  visNetwork             2.0.5     2018-12-05 [1] CRAN (R 3.4.4)                          
##  vsn                  * 3.46.0    2019-01-04 [1] Bioconductor                            
##  withr                  2.1.2     2018-03-15 [1] CRAN (R 3.5.1)                          
##  xfun                   0.4       2018-10-23 [1] CRAN (R 3.4.4)                          
##  XML                    3.98-1.16 2018-08-19 [1] CRAN (R 3.4.4)                          
##  xml2                   1.2.0     2018-01-24 [1] CRAN (R 3.4.4)                          
##  xtable                 1.8-3     2018-08-29 [1] CRAN (R 3.4.4)                          
##  XVector                0.20.0    2018-09-25 [1] Bioconductor                            
##  yaml                   2.2.0     2018-07-25 [1] CRAN (R 3.5.1)                          
##  zip                    1.0.0     2017-04-25 [1] CRAN (R 3.5.1)                          
##  zlibbioc               1.26.0    2018-09-25 [1] Bioconductor                            
##  zoo                    1.8-4     2018-09-19 [1] CRAN (R 3.4.4)                          
## 
## [1] /home/angel/R/x86_64-pc-linux-gnu-library/3.5
## [2] /usr/local/lib/R/site-library
## [3] /usr/lib/R/site-library
## [4] /usr/lib/R/library
```

## References
