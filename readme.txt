README: CVD Burden Decomposition Script

OVERVIEW
This R script accompanies the following manuscript:
Global Burden of Cardiovascular Diseases and Risks 2023 Collaborators. Global, regional, and national burden of cardiovascular diseases and risk factors in 204 countries and territories, 1990-2023. Journal of the American College of Cardiology. 24 September 2025.

Both the manuscript and this code are compliant with the Guidelines for Accurate and Transparent Health Estimates Reporting (GATHER), ensuring clear documentation of data sources, analytic methods, and reproducibility of results.

The script performs a decomposition analysis of changes in disease burden (e.g., DALYs, deaths, YLDs) for cardiovascular diseases (CVD) over time, attributing changes to four main factors:
1. Population growth
2. Population aging (age structure)
3. Risk-deleted rate changes
4. Risk exposure changes

The decomposition is performed for a specified location and measure, using uncertainty propagation via draws. The script is designed to run in a high-performance computing (HPC) environment and is modular for use with the Global Burden of Disease (GBD) data pipeline. The code to produce CVD estimates (death, DALYs, prevalence) by age, sex, year, and location will be made publicly available upon publication of the GBD 2023 enterprise. Please see the IHME website for further information: https://www.healthdata.org/research-analysis/gbd


MAIN STEPS
1. Setup: Loads required libraries, defines custom operators, and sources utility functions.
2. Initialize variables: Sets global variables for location, measure, years, release versions, and IDs.
3. Data loading & preparation: Loads cause, risk, population, and mediation metadata; pulls draws of risk-attributable burden, population, and rates.
4. Decomposition analysis: Calculates the contribution of each factor to changes in disease burden using Das Gupta’s decomposition method.
5. Mediation adjustment: Applies mediation matrix to allocate overlapping risk effects.
6. Raking & scaling: Ensures decomposed factors sum to total observed change, raking as needed.
7. Uncertainty quantification: Propagates uncertainty via draws, calculates lower/upper UI bounds.
8. Output formatting: Aggregates and saves results including both percent and count outputs.


INPUTS
1. GBD functions: Functions like get_draws, get_population, get_demographics, etc., are assumed to be in the sourced utility scripts.
2. Metadata csv's: Mediation matrix (mediation_matrix_draw_gbd_2023.csv), strategy flags (strategies_for_pafsofone.csv), and REI metadata (cfr_sev.csv) are loaded as needed.
3. Directory paths: Paths are set via global variables (PATH, outdir) and must be edited for your environment.
4. Global variables: Set for location ID, measure ID, years, release and run IDs.


OUTPUTS
1. Decomposition results: Excel files containing decomposition results with percent and count change by factor, CVD cause, risk, and uncertainty intervals.
2. Output is saved to: measure_<MEASURE_ID>/<LOCATION_ID>_decomp.xlsx in the specified output directory.

Example output columns:
cause_id, rei_id, sex_id, factor_format ("population_effect", "age_structure_effect", etc.)
change (mean percent change), lower_change, upper_change
change_count (mean count change), lower_change_count, upper_change_count
total (total percent change), lower_total, upper_total


USAGE
1. Configure file paths: Edit all PATH and file locations to point to your actual directories and metadata files.
2. Set global variables: Update IDs for location, measure, release, etc. 
3. Source utility functions: Place all required utility functions in the specified central directory.
4. Run script: Execute the script in R. For HPC, submit as a SLURM array job using the SLURM_ARRAY_TASK_ID environment variable.


DEPENDENCIES
1. R version 3.6 or more recent recommended
2. CRAN packages: data.table, dplyr, plyr, zoo, magrittr, parallel, rhdf5, openxlsx
3. IHME GBD utility functions (must be sourced from your central repo)
4. Metadata CSV files in specified locations


KEY FUNCTIONS
1. loadData(): Loads and merges GBD draws for risks, causes, population, and rates; prepares main data table.
2. decompFactors(): Implements Das Gupta's decomposition to split total change into population, age, rate, and risk effects.
3. meatiator(): Applies mediation matrix to adjust overlapping risk effects.
4. riskRaker(): Rakes and scales decomposed factors to ensure consistency with total observed change.
5. Uncertainty calculation: Quantiles calculated for lower andupper bounds (2.5% and 97.5%) across draws.
6. Output formatting: Aggregates by cause/risk, merges with metadata, and saves as Excel.