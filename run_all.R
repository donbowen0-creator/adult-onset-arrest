###############################################################################
# run_all.R  --  RUN THIS ONE FILE
#
# Rebuilds the whole analysis from the raw ICPSR files:
#   00_setup.R             packages, folders, and every up-front decision
#   01_build_data.R        builds variables; prints N after every filter
#   02_models.R            self-reported onset models (companion analysis)
#   03_tables.R            tables for the self-report analysis
#   04_figures.R           figures for the self-report analysis
#   05_arrest_onset.R      adult-onset ARREST: two stacked risk sets
#   05b_revision_checks.R  own substance use, significance tests, no-imputed check
#   05c_sex_interactions.R sex interactions and sex-stratified models
#   06_manuscript_tables.R journal tables and figure files
# then saves the R and package versions (sessionInfo) for the Methods section.
#
# HOW TO RUN (RStudio):
#   1. Open this file in RStudio.
#   2. Click "Source" (top right of the script pane).
# This file now finds its own folder, so you no longer need to set the
# working directory by hand. If it cannot work out where it is, it says so
# and tells you what to do.
###############################################################################

start_time <- Sys.time()

## ---------------------------------------------------------------------------
## Find the folder this file is in, and work from there
## ---------------------------------------------------------------------------
find_script_dir <- function() {
  # 1. when the file is sourced, R records its path
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  # 2. when it is run with Rscript
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  # 3. when it is run from the RStudio editor
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  NULL
}

script_dir <- find_script_dir()
if (!is.null(script_dir)) setwd(script_dir)

stages <- c("00_setup.R", "01_build_data.R", "02_models.R", "03_tables.R",
            "04_figures.R", "05_arrest_onset.R", "05b_revision_checks.R",
            "05c_sex_interactions.R", "06_manuscript_tables.R")
missing <- stages[!file.exists(stages)]
if (length(missing)) {
  stop("\n\nR is looking for the scripts in:\n   ", getwd(),
       "\nbut cannot find:\n   ", paste(missing, collapse = "\n   "),
       "\n\nFix it in one of these ways:",
       "\n  * In RStudio: Session > Set Working Directory > To Source File Location,",
       "\n    then click Source again.",
       "\n  * Or edit the line below and run it, using forward slashes:",
       "\n      setwd(\"C:/path/to/AddHealth_Onset_Rebuild\")",
       "\nFiles R can see in the current folder:\n   ",
       paste(head(list.files(), 25), collapse = "\n   "), "\n", call. = FALSE)
}

## ---------------------------------------------------------------------------
## Run every stage, logging to screen and to file
## ---------------------------------------------------------------------------
dir.create("output/logs", showWarnings = FALSE, recursive = TRUE)
log_con <- file("output/logs/run_log.txt", open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("Add Health rebuild -- run started", format(start_time), "\n")
cat("Working folder:", getwd(), "\n")

ok <- TRUE
for (s in stages) {
  cat("\n>>>", s, "\n")
  ok <- tryCatch({ source(s, echo = FALSE); TRUE },
    error = function(e) {
      cat("\n\n*** THE RUN STOPPED IN", s, "***\n", conditionMessage(e),
          "\nThe full log is in output/logs/run_log.txt -- send that file along",
          "when asking for help.\n")
      FALSE
    })
  if (!ok) break
}

writeLines(capture.output(sessionInfo()), "output/logs/sessionInfo.txt")
cat("\nR and package versions saved to output/logs/sessionInfo.txt\n")
cat("Run finished", format(Sys.time()), "-- elapsed",
    round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "min\n")
if (ok) {
  cat("\nALL STAGES COMPLETED.\n",
      "  Self-report analysis tables: output/tables/ALL_TABLES.html\n",
      "  Arrest analysis tables:      output/arrest_onset/ARREST_ONSET.html\n",
      "  Manuscript tables and figures: output/manuscript/\n")
}

sink(type = "message"); sink(); close(log_con)
