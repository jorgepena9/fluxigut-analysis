# Bidirectional interactions between the gut microbiota and fluorochemical biotransformation and bioactivity

This study investigated the effects of 15 structurally diverse fluorinated chemicals on the composition of the human gut microbiome using an ex vivo high-throughput cultivation system. Fecal samples from 13 healthy human donors were incubated with the different compounds for 48 hours. 16S rRNA gene amplicon sequencing (V4 region, 515F/806R primers) was performed on an Illumina NextSeq 2000 using paired-end reads. This dataset demonstrates that exposure to fluorinated chemicals can induce compound-specific shifts in microbial diversity and community composition, highlighting their capacity to alter gut microbial ecology.

---

## Repository structure

```
fluxigut-analysis/
├── README.md
├── .gitignore
├── Fluxigut_analysis.Rmd                # Main R analysis 
├── Fluxigut_qiime2_processing.ipynb     # QIIME2 preprocessing pipeline
├── data/
│   ├── metadata/
│   │   └── Fluxigut_metadata_clean.csv  # Sample metadata (13 donors, 16 treatments)
│   ├── sequencing/
│   │   ├── feature-table.tsv            # ASV feature table (QIIME2 output)
│   │   └── taxonomy_SILVA_weighted.tsv  # SILVA 138 weighted taxonomy
│   ├── concentrations/
│   │   └── RP_DataClean.csv             # Targeted metabolomics (drug concentrations)
│   └── od600/
│       └── OD_summary.csv               # OD600 data
└── results/                             # Figures 
```

---

### Raw sequencing data

Raw 16S rRNA amplicon sequencing data (V4 region, 515F/806R, Illumina NextSeq 2000) are deposited at the European Nucleotide Archive (ENA) under accession **PRJEB112299**.

### QIIME2 preprocessing

The notebook `code/Fluxigut_qiime2_processing.ipynb` documents the full preprocessing pipeline from raw reads to the feature table and taxonomy files in `data/sequencing/`. 

### R analysis

The Rmd `code/Fluxigut_analysis.Rmd` runs the full downstream analysis and produces all figures in the manuscript. 

