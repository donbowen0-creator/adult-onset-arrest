# Adult-onset offending and arrest in Add Health: analysis code

R code that builds the analysis files and estimates every model in a study of adult-onset offending and first arrest, using the public-use data of the National Longitudinal Study of Adolescent to Adult Health (Add Health).

**The data are not included here.** Add Health public-use data are distributed by ICPSR under terms that do not permit redistribution. To run this code you need the ICPSR 21600 public-use download, which is free once you register.

## What the code does

Two analyses run from the same build.

**Self-reported onset (stages 1 to 4).** Adolescent delinquency at Wave I, then adult-onset offending and illicit drug use at Waves III, IV and V, among respondents who reported no delinquency and no illicit drug use as adolescents. Weighted logistic models with cluster-robust standard errors.

**First arrest (stages 5, 5b, 5c and 6).** Age at first self-reported arrest among respondents with no arrest before 18, modeled as a discrete-time hazard across two stacked risk sets divided at the Wave III interview, so that adult role measures always precede the outcome. Includes the composition of the adult-onset category under two screens of juvenile offending, tests separating detection from involvement, a split of the peer measure into substance use and co-offending, and tests of sex differences.

Every analytic decision, including which items count as offending in each wave, which weight goes with which model, and how the risk sets are defined, is set in one place: section 3 of `00_setup.R`. Change a decision there and the whole pipeline follows it.

## Getting set up

1. Install R from https://cran.r-project.org and, if you want a friendlier interface, RStudio Desktop.
2. Download the Add Health public-use data (ICPSR 21600) from ICPSR. Keep the zip file; there is no need to unzip it.
3. Open `00_setup.R` and set `icpsr_source` near the top to the location of that zip file on your computer, using forward slashes, for example `C:/Users/you/Desktop/ICPSR_21600-V26.zip`. As a fallback the scripts will also look for the zip in a `data_raw` folder inside the project.
4. Open `run_all.R` and click Source. It finds its own folder, installs the three packages it needs (haven, survey, ggplot2) the first time, extracts the handful of ICPSR files the analysis uses, and runs everything in order.

A full run takes two to three minutes on an ordinary laptop and is designed to work within 4 GB of memory.

## What you get

```
output/
  logs/        run log, sessionInfo, sample flow at every filtering step
  tables/      self-report analysis tables (CSV, plus one HTML file)
  figures/     self-report analysis figures (PNG and PDF)
  arrest_onset/  first-arrest tables, robustness checks, figures (CSV, HTML, PNG)
  manuscript/  journal-ready tables and figure files
```

The run prints the sample size after every filtering step and, at the end of stage 5, prints the headline results, so the log alone records what a run produced.

## The scripts

| File | What it does |
| --- | --- |
| `00_setup.R` | Packages, folders, file locations, and every analytic decision |
| `01_build_data.R` | Builds all variables from Waves I, III, IV and V; prints the sample flow |
| `02_models.R` | Self-reported onset models, weighted and unweighted |
| `03_tables.R` | Tables for the self-report analysis |
| `04_figures.R` | Figures for the self-report analysis |
| `05_arrest_onset.R` | First-arrest risk sets, discrete-time hazard models, composition of the category |
| `05b_revision_checks.R` | Own substance use as a control, prevalence-ratio tests, imputed-age check |
| `05c_sex_interactions.R` | Sex interactions and sex-stratified models |
| `06_manuscript_tables.R` | Journal-ready tables and figure files |
| `run_all.R` | Runs every stage in order and saves the R and package versions |

## Requirements

R 4.3 or later, with haven, survey and ggplot2. The scripts install them on first run if they are missing.

## Data citation

Hummer, R. A., Aiello, A. E., Harris, K. M., & Udry, J. R. *National Longitudinal Study of Adolescent to Adult Health (Add Health), 1994-2025 [Public Use]* (Version V26). Inter-university Consortium for Political and Social Research. https://doi.org/10.3886/ICPSR21600.v26

This research uses data from Add Health, funded by grant P01 HD31921 from the Eunice Kennedy Shriver National Institute of Child Health and Human Development, with cooperative funding from 23 other federal agencies and foundations. No direct support was received from that grant for this analysis.

## License

MIT. See `LICENSE`.

## Contact

Donald Bowen, DBA
don@cemresearchlab.com
https://www.cemresearchlab.com
ORCID: https://orcid.org/0000-0003-3527-9462
