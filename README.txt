ADD HEALTH ADULT-ONSET ANALYSIS -- REBUILT FROM SCRATCH IN R
=============================================================

WHAT IS IN THIS FOLDER
  run_all.R          <- the only file you need to run
  00_setup.R         packages, folders, and EVERY analytic decision (read this)
  01_build_data.R    Stage 1: builds variables, prints N after every filter
  02_models.R        Stage 2: weighted + unweighted models, strain x peer
  03_tables.R        Stage 3: all tables (CSV + one HTML file)
  04_figures.R       Stage 4: all figures (PNG + PDF)
  05_arrest_onset.R  Stage 5: adult-onset ARREST (age at first arrest) as
                     two stacked risk sets -- window 1 runs from age 18 to
                     the Wave III interview, window 2 from Wave III to the
                     last interview with lagged Wave III roles. Discrete-time
                     survival on person-years. Results in output/arrest_onset/
                     (open ARREST_ONSET.html).
  06_manuscript_tables.R
                     Stage 6: journal-ready tables (CSV) and figure files
                     (Figure_1, Figure_2 as PDF + 600-dpi PNG) for the
                     detection paper, in output/manuscript/.
  manuscript_build/  make_tables_docx.js turns the Stage 6 CSVs into
                     Tables.docx (run: node make_tables_docx.js
                     ../output/manuscript Tables.docx). Needs Node.js; optional.
  data_raw/          put ICPSR_21600-V26.zip here (you already have it)
  output/            results from the test run (they are regenerated each run)

ONE-TIME SETUP
  1. Install R (https://cran.r-project.org) and RStudio Desktop (free).
  2. The scripts are set to read the zip straight from your Desktop:
        C:/Users/don/OneDrive/Desktop/ICPSR_21600-V26.zip
     If you move it, change icpsr_source at the top of 00_setup.R (use
     forward slashes). As a fallback you can also copy the zip into the
     project's data_raw folder. Do NOT unzip it --
     the script pulls out only the files it needs (keeps memory low
     on a 4 GB laptop). An already-unzipped ICPSR_21600 folder also works.

EACH TIME YOU RUN IT
  1. Open run_all.R in RStudio.
  2. Menu: Session > Set Working Directory > To Source File Location
  3. Click "Source" (top right of the script window).
  4. Wait for "ALL STAGES COMPLETED" (about 1-3 minutes; the first run
     also installs three packages and unzips the data).
  5. Open output/tables/ALL_TABLES.html (double-click) to see every table.

  If it stops with an error, the complete record is in
  output/logs/run_log.txt -- send that file along when asking for help.

WHERE THINGS LAND
  output/logs/run_log.txt        everything printed during the run,
                                 including the sample size at every step
  output/logs/sample_flow.csv    the N-by-step table (Table S1)
  output/logs/sessionInfo.txt    R and package versions for Methods
  output/tables/                 Tables 1, 1b, 2a, 2b, 3, 4, S1, S2
  output/figures/                Figures 1-5

THE DECISIONS (full detail and item lists are in 00_setup.R, section 3)
  Waves used: I, III, IV, V. Wave II is not used. Arrest is not used anywhere.
  Wave I delinquency: any of 13 criminal acts in the past 12 months (lying
    to parents, running away, and loud/rowdy excluded).
  Adult onset = self-reported offending OR illicit drug use (marijuana,
    cocaine, meth, heroin, other illegal drugs, nonmedical prescription
    drugs at Wave V). Alcohol is not counted.
  Risk set: no Wave I delinquency, no Wave I illicit drug use, and no first
    drug use before 18 (Wave IV age-at-first-use).
  PRIMARY, all adult onset: Wave III (age 19+) is an onset window, so onset
    can occur at Waves III, IV or V (~ages 19-43). Weight GSW1345.
  SECONDARY, late adult onset: Wave III is a screen instead (no offending or
    drug use by Wave III); onset at Waves IV-V (~ages 25-43). Weight GSW1345.
  SENSITIVITY: Wave I screen only. Weight GSW145.
  Offending onset and drug-use onset are also modeled separately, plus
    drug-use onset excluding marijuana.
  Models: A core; B + strain x peer; C + physical abuse, sexual abuse
    (before 18, Wave IV), childhood ADHD symptoms (Wave III recall);
    D + Wave I binge drinking (robustness).
  Weights: GSWGT1 for H1; clustering on CLUSTER2 throughout. Weighted is
    primary; unweighted keeps the cluster correction.

TO CHANGE A DECISION
  Edit the item list or switch in 00_setup.R, section 3, save, and run
  run_all.R again. The sample-flow table will show exactly how the Ns move.
