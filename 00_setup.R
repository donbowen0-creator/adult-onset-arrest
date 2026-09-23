###############################################################################
# 00_setup.R
# Add Health adult-onset analysis -- SETUP AND UP-FRONT DECISIONS
#
# This file does three things:
#   1. Loads (and if needed installs) the R packages the project uses.
#   2. Sets up folder paths and finds the ICPSR 21600 data.
#   3. Writes down, in one place, EVERY analytic decision:
#        - which items count as delinquency at Wave I
#        - which items count as offending / arrest at Waves IV and V
#        - which weight goes with which model
#      Change a decision HERE and every later script follows it.
#
# You never need to run this file by itself. run_all.R runs it first.
###############################################################################

# ---------------------------------------------------------------------------
# 1. PACKAGES
# ---------------------------------------------------------------------------
# haven   = reads the Stata (.dta) and SAS (.sas7bdat) files ICPSR supplies
# survey  = survey-weighted models with clustered standard errors
# ggplot2 = figures
needed <- c("haven", "survey", "ggplot2")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing package: ", p, " (one-time, needs internet)")
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(ggplot2)
})

# Clusters with only one case in a subpopulation should not stop the models.
# "adjust" centres those clusters on the grand mean (a conservative, standard
# choice). This affects standard errors only, never coefficients.
options(survey.lonely.psu = "adjust")

# ---------------------------------------------------------------------------
# 2. FOLDERS AND DATA LOCATION
# ---------------------------------------------------------------------------
# The project folder is wherever these scripts live. run_all.R sets the
# working directory to that folder before sourcing this file.
proj_dir <- getwd()
raw_dir  <- file.path(proj_dir, "data_raw")      # put the ICPSR zip here
work_dir <- file.path(proj_dir, "data_work")     # built data files go here
out_dir  <- file.path(proj_dir, "output")
tab_dir  <- file.path(out_dir, "tables")
fig_dir  <- file.path(out_dir, "figures")
log_dir  <- file.path(out_dir, "logs")
for (d in c(raw_dir, work_dir, out_dir, tab_dir, fig_dir, log_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# WHERE YOUR ICPSR DOWNLOAD LIVES ------------------------------------------
# Point this at the ICPSR zip (or at an already-unzipped ICPSR_21600 folder).
# Use forward slashes, even on Windows -- R reads "C:/Users/..." fine but
# treats a single backslash as an escape character.
# If this path does not exist, the script falls back to any ICPSR_21600*.zip
# sitting in the project's data_raw folder.
icpsr_source <- "C:/Users/don/OneDrive/Desktop/ICPSR_21600-V26.zip"

# The ICPSR files this project needs (and nothing else, to save memory).
icpsr_files <- c(
  w1_main    = "ICPSR_21600/DS0001/21600-0001-Data.dta",   # Wave I in-home
  w1_weight  = "ICPSR_21600/DS0004/21600-0004-Data.dta",   # Wave I weights
  w1_net     = "ICPSR_21600/DS0003/21600-0003-Data.dta",   # Wave I network variables
  w3_main    = "ICPSR_21600/DS0008/21600-0008-Data.dta",   # Wave III in-home
  w4_main    = "ICPSR_21600/DS0022/21600-0022-Data.dta",   # Wave IV in-home
  w4_weight  = "ICPSR_21600/DS0031/21600-0031-Data.dta",   # Wave IV weights
  w5_main    = "ICPSR_21600/DS0032/21600-0032-Data.dta",   # Wave V survey
  w5_zip     = "ICPSR_21600/DS0042/21600-0042-Zipped_package-MULTI.zip"
)

# Finds the data whether you unzipped the ICPSR download or not.
# Option A: data_raw/ICPSR_21600/...   (already unzipped)
# Option B: data_raw/ICPSR_21600-V26.zip (or any ICPSR_21600*.zip) -- the
#           script pulls out only the files it needs, one time.
locate_icpsr <- function() {
  have_all <- all(file.exists(file.path(raw_dir, icpsr_files)))
  if (!have_all) {
    zips <- list.files(raw_dir, pattern = "^ICPSR_21600.*\\.zip$",
                       full.names = TRUE)
    if (length(zips) == 0) {
      stop("\n\nCannot find the Add Health data.\n",
           "Put ICPSR_21600-V26.zip (or the unzipped ICPSR_21600 folder)\n",
           "inside this folder:\n   ", raw_dir, "\n")
    }
    message("Extracting the needed files from ", basename(zips[1]),
            " (one time only; takes a minute)...")
    unzip(zips[1], files = icpsr_files, exdir = raw_dir)
  }
  # Wave V weights are inside a second zip inside DS0042.
  w5w <- file.path(raw_dir, "ICPSR_21600/DS0042/p5weight.sas7bdat")
  if (!file.exists(w5w)) {
    unzip(file.path(raw_dir, icpsr_files[["w5_zip"]]),
          files = "p5weight.sas7bdat",
          exdir = file.path(raw_dir, "ICPSR_21600/DS0042"))
  }
  paths <- as.list(file.path(raw_dir, icpsr_files))
  names(paths) <- names(icpsr_files)
  paths$w5_weight <- w5w
  paths
}

# ---------------------------------------------------------------------------
# 3. UP-FRONT DECISIONS  (the part to read carefully and cite in Methods)
# ---------------------------------------------------------------------------

## 3a. WAVE I DELINQUENCY (the H1 outcome, and the screen for "abstainers")
# Past-12-month self-report items H1DS1-H1DS15 (0 = never ... 3 = 5+ times).
# INCLUDED: the 12 items describing acts that would be crimes at any age.
# EXCLUDED: H1DS3 (lied to parents), H1DS7 (ran away), H1DS15 (loud/rowdy in
#           public) -- status or nuisance behaviours, not offending.
# PLUS:     H1FV7 (pulled a knife/gun on someone), so that every act counted
#           as adult offending at Wave V (which includes that act) was also
#           screened for at Wave I.
w1_delinq_items <- c(
  H1DS1  = "painted graffiti",
  H1DS2  = "damaged property",
  H1DS4  = "shoplifted",
  H1DS5  = "serious physical fight",
  H1DS6  = "hurt someone badly enough to need care",
  H1DS8  = "took a car without permission",
  H1DS9  = "stole something worth > $50",
  H1DS10 = "burglary (went into house/building to steal)",
  H1DS11 = "used/threatened weapon to get something",
  H1DS12 = "sold drugs",
  H1DS13 = "stole something worth < $50",
  H1DS14 = "group fight",
  H1FV7  = "pulled a knife or gun on someone"
)

## 3b. WAVE IV SELF-REPORTED OFFENDING (past 12 months, ages ~24-32)
# INCLUDED: acts that parallel the Wave I list.
# EXCLUDED: H4DS9 (used someone's credit card) and H4DS10 (bad check) --
#           fraud has no Wave I counterpart, so including it would create
#           "onset" partly by measurement. H4DS13-H4DS18 are victimization
#           items, not offending.
w4_offend_items <- c(
  H4DS1  = "damaged property",
  H4DS2  = "stole something worth > $50",
  H4DS3  = "burglary",
  H4DS4  = "used/threatened weapon to get something",
  H4DS5  = "sold drugs",
  H4DS6  = "stole something worth < $50",
  H4DS7  = "group fight",
  H4DS8  = "bought/sold/held stolen property",
  H4DS11 = "serious physical fight",
  H4DS12 = "hurt someone badly (asked only if fought; skip = 0)",
  H4DS19 = "pulled a knife or gun on someone",
  H4DS20 = "shot or stabbed someone"
)

## 3c. WAVE V SELF-REPORTED OFFENDING (past 12 months, ages ~33-43)
# Wave V asked only these five offending items; all five are used.
w5_offend_items <- c(
  H5CJ1A = "damaged property",
  H5CJ1B = "stole something worth > $50",
  H5CJ1C = "sold drugs",
  H5CJ1D = "physical fight",
  H5CJ1E = "pulled a knife or gun on someone"
)

## 3d. ARREST -- NOT USED
# Arrest plays no role in any outcome or screen. Onset is defined by the
# respondent's own reported behaviour: offending and illicit drug use.

## 3e. MISSING-DATA RULE FOR MULTI-ITEM MEASURES
# A "none of these acts" (0) score requires at least this share of the items
# to have a real answer. Any single "yes" makes the score 1 regardless.
min_valid_share <- 0.80

## 3f. WAVE III (2001-02, ages ~18-26): offending
# Self-reported offending, past 12 months: acts parallel to Waves I and IV.
# (Fraud, gun-carrying, gang, and handgun-ownership items excluded.)
w3_offend_items <- c(
  H3DS1 = "damaged property",       H3DS2 = "stole something worth > $50",
  H3DS3 = "burglary",               H3DS4 = "used/threatened weapon to get something",
  H3DS5 = "sold drugs",             H3DS6 = "stole something worth < $50",
  H3DS7 = "group fight",            H3DS8 = "bought/sold/held stolen property",
  H3DS11 = "used weapon in a fight", H3DS17 = "hurt someone badly (count)",
  H3DS18H = "pulled a knife or gun on someone", H3DS18I = "shot or stabbed someone"
)
# Wave III acts count as ADULT only if the respondent was 19+ at Wave III
# (so the whole 12-month window falls at 18+).

## 3g. ILLICIT DRUG USE AS CRIMINAL ACTIVITY ("substance onset")
# Illicit drug use is itself an offense, so it counts as adult onset.
# Alcohol is NOT included (legal for adults). Measures, by wave:
#   Wave I  (screen): ever used marijuana, cocaine, inhalants, or other
#           illegal drugs (H1TO30, H1TO34, H1TO37, H1TO40). Users are not
#           abstainers and are removed from the risk set.
#   Wave III: past-year marijuana, cocaine, crystal meth, other drugs
#           (H3TO109/112/115/118); used as an onset window only at age 19+.
#           Any use since June 1995 incl. nonmedical prescription drugs
#           (H3TO105A-D, H3TO108/111/114/117/120) is the late-onset screen.
#   Wave IV: past-year marijuana (H4TO70) or past-year use of the most-used
#           other illegal drug (H4TO98). Age at first use (H4TO68, H4TO96)
#           removes anyone who FIRST used before 18.
#           (Wave IV nonmedical prescription use is "ever", with no age, so
#           its timing cannot be placed; it is not used.)
#   Wave V:  past-30-day marijuana (H5TO21), first-ever marijuana use between
#           Waves IV and V (H4TO65B = no, H5TO20 = yes), past-30-day cocaine,
#           meth, heroin, other drugs (H5TO27A-D), and nonmedical
#           prescription drug use (H5TO26A-D).
# SENSITIVITY: the same outcome EXCLUDING marijuana, because recreational
# marijuana was legal for adults in several states during Wave V.

## 3h. ADULT-ONSET OUTCOMES
# Risk set: no Wave I delinquency, no Wave I illicit drug use, and no first
# drug use before 18 reported at Wave IV.
# Any adult onset = self-reported offending OR illicit drug use.
# The two components are also modeled separately.
# PRIMARY  -- all adult onset: Wave III is an onset WINDOW (age 19+), so
#     onset can occur at Waves III, IV or V (~ages 19-43).
# SECONDARY -- late adult onset: Wave III is a SCREEN instead (no offending
#     and no drug use by Wave III); onset at Waves IV-V (~ages 25-43). This
#     is the stricter design, but it leaves few self-reported offending
#     cases, so it is reported as a robustness check.
# SENSITIVITY -- Wave I screen only (no Wave III); onset at Waves IV-V.
# Wave II is NOT used.

## 3i. WEIGHTS -- decided before any model is run
# Model family                     Waves used        Weight    Cluster
# H1  (Wave I delinquency)         I                 GSWGT1    CLUSTER2
# H2  primary (late onset)         I, III, IV, V     GSW1345   CLUSTER2
# H2  secondary (all adult onset)  I, III, IV, V     GSW1345   CLUSTER2
# H2  sensitivity (W1 screen)      I, IV, V          GSW145    CLUSTER2
# The public-use files carry no stratum variable, so designs are
# cluster-only (ids = ~CLUSTER2), the standard public-use specification.
# Weighted results are PRIMARY. Unweighted results keep the cluster
# correction so that the ONLY difference between the two is the weights.
weight_plan <- data.frame(
  model   = c("H1: Wave I delinquency",
              "H2 primary: all adult onset (W3 window; onset W3-W5)",
              "H2 secondary: late adult onset (W3 screen; onset W4-W5)",
              "H2 sensitivity: earlier design (W1 screen only)"),
  waves   = c("I", "I, III, IV, V", "I, III, IV, V", "I, IV, V"),
  weight  = c("GSWGT1", "GSW1345", "GSW1345", "GSW145"),
  cluster = "CLUSTER2",
  primary = "weighted"
)

## 3j. PREDICTORS (all measured at Wave I)
# Strain       = CES-D, 19 items (H1FS1-H1FS19), 0-3 each; H1FS4, H1FS8,
#                H1FS11, H1FS15 reverse-coded. Mean of answered items x 19,
#                if >= 15 items answered.
# Peer deviance= number of 3 best friends who smoke (H1TO9), drink monthly
#                (H1TO29), use marijuana monthly (H1TO33); 0-9, all 3 needed.
# Bond deficit = mean of z-scored items, coded so HIGHER = WEAKER bonds:
#                school: close to people (H1ED19), part of school (H1ED20),
#                happy at school (H1ED22); family: parents care (H1PR3),
#                family understands (H1PR5), family has fun (H1PR7),
#                family pays attention (H1PR8). Needs >= 4 of 7 items
#                (so youths not in school are scored on the family items).
# Strain, peer deviance, and bond deficit enter models as z-scores
# (standardized on everyone with a valid Wave I score), so the strain x peer
# interaction is interpretable and comparable across models.
cesd_items   <- paste0("H1FS", 1:19)
cesd_reverse <- c("H1FS4", "H1FS8", "H1FS11", "H1FS15")
peer_items   <- c("H1TO9", "H1TO29", "H1TO33")
school_items <- c("H1ED19", "H1ED20", "H1ED22")      # 1 = strongly agree
family_items <- c("H1PR3", "H1PR5", "H1PR7", "H1PR8") # 5 = very much

## 3k. CONTROLS (Wave I)
# age at interview, male, race/ethnicity (Hispanic; non-Hispanic Black;
# non-Hispanic other; reference = non-Hispanic White), and any resident
# parent with a 4-year college degree or more. About 5% of respondents could
# not report a resident parent's education (don't know / no resident parent);
# rather than drop them, they are scored 0 on par_college and flagged with
# par_ed_unknown = 1 (a missing-indicator), so the sample is not cut by ~350.
controls <- c("age_w1", "male", "hisp", "black", "other_race",
              "par_college", "par_ed_unknown")

## 3l. EARLY ADVERSITY BLOCK (added in Model C of the onset models)
# Physical abuse before 18 (H4MA3): a parent or adult caregiver hit, kicked,
#   or threw you at least once (1 = yes). Reported at Wave IV.
# Sexual abuse before 18 (H4MA5): a parent or adult caregiver touched you
#   sexually or forced sexual contact at least once (1 = yes). Wave IV.
# Childhood ADHD symptoms (H3RA1-H3RA17), recalled at Wave III for ages
#   5-12: number of the 17 DSM-IV symptoms rated "often" or "very often"
#   (needs >= 14 answered; prorated to 17). Enters models as a z-score.
#   "Probable ADHD" (>= 6 inattentive OR >= 6 hyperactive-impulsive
#   symptoms) is reported in the descriptives.
# All three are retrospective; they are measured before the onset window
# but reported during adulthood.
adhd_inatt <- paste0("H3RA", c(1, 3, 5, 7, 9, 11, 13, 15, 17))
adhd_hyper <- paste0("H3RA", c(2, 4, 6, 8, 10, 12, 14, 16))
adversity  <- c("phys_abuse", "sex_abuse", "adhd")

## 3m. ROBUSTNESS: OWN WAVE I BINGE DRINKING
# Everyone in the onset samples is a Wave I non-user of illicit drugs, but
# some drank. Because the peer-deviance scale is built from FRIENDS' substance
# use, a robustness model (Model D) adds own past-year binge drinking at
# Wave I (H1TO17; 1 = any occasion of 5+ drinks).

# ---------------------------------------------------------------------------
# 4. SMALL HELPER FUNCTIONS used by the other scripts
# ---------------------------------------------------------------------------

# Turn an ICPSR labelled column into plain numbers and blank out codes
# that mean refused / don't know / not applicable.
clean_num <- function(x, missing_codes) {
  x <- as.numeric(haven::zap_labels(x))
  x[x %in% missing_codes] <- NA
  x
}

# Any-act indicator: 1 if any item > 0; 0 if no item > 0 AND enough items
# were answered; otherwise NA.
any_act <- function(df, items, min_share = min_valid_share) {
  m <- as.matrix(df[, items, drop = FALSE])
  n_valid <- rowSums(!is.na(m))
  any_yes <- rowSums(m > 0, na.rm = TRUE) > 0
  out <- ifelse(any_yes, 1,
                ifelse(n_valid >= ceiling(min_share * length(items)), 0, NA))
  as.numeric(out)
}

# Sample-flow recorder: every filter step is printed AND saved.
flow_log <- data.frame(sample = character(), step = character(),
                       n = integer(), dropped = integer())
record_step <- function(sample, step, n) {
  prev <- flow_log$n[flow_log$sample == sample]
  dropped <- if (length(prev)) tail(prev, 1) - n else NA
  flow_log[nrow(flow_log) + 1, ] <<- list(sample, step, n, dropped)
  cat(sprintf("  [%s] %-62s N = %5d%s\n", sample, step, n,
              if (is.na(dropped)) "" else sprintf("  (dropped %d)", dropped)))
}

# z-score helper
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

cat("Setup complete. Project folder:", proj_dir, "\n")
