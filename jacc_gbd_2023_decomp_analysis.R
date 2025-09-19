rm(list=ls())

# Set up environment -----------------------------------------------------------
library(data.table)
library(plyr)
library(zoo)
library(magrittr)
library(parallel)
library(rhdf5)
library(dplyr)
'%ni%' <- Negate('%in%')

# Define functions -------------------------------------------------------------
central <- "PATH"
for (func in paste0(central, list.files(central))) source(paste0(func))

# Define global variables ------------------------------------------------------
args <- commandArgs(trailingOnly = T)

## Retrieving array task_id
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(task_id)

location_id <- "ID" ## global
measure_id <- "ID"  ## DALYs

print(paste0("Creating decomp analysis for  ", location_id, ' and measure_id ', measure_id))
year_ids <- c(1990, 2023) 
release_id <- "ID"
ndraws <- 250 

cvid <- "ID" 
burdenator_id <- "ID"  
codcorrect_id <- 'ID' 
dalynator_id <- 'ID'   
sev_cvid <- 'ID'     
population_estimate_version <- 'ID'

outdir <- 'PATH'

demo <- get_demographics(gbd_team="epi", release_id = release_id)
age_group_ids <- demo$age_group_id
sex_ids <- demo$sex_id
reimetadf <- get_rei_metadata(rei_set_id = 1, release_id = release_id)
reis <- copy(reimetadf)

# set up REI IDs
# for lbw/sg and air we estimate the joint directly
reis <- reis[!rei_id %in% c(86, 87), ]
reis[rei_id %in% c(380), most_detailed := 1]
parent_cols <- paste0("parent_",min(reis$level):max(reis$level))
reis <- reis[, .(rei_id, path_to_top_parent, most_detailed)]
reis[, (parent_cols) := tstrsplit(path_to_top_parent, ",")][, path_to_top_parent := NULL]
reis <- data.table::melt(reis, id.vars = c("rei_id", "most_detailed"), value.vars = parent_cols,
                         value.name = "parent_id", variable.name = "level")
reis[, parent_id := as.integer(parent_id)]
reis <- reis[!is.na(parent_id) & (most_detailed == 1 | rei_id == parent_id), ]
reis[, level := NULL]
cause_ids <- get_ids('cause')
cvd_cause_ids <- cause_ids[grep('cvd',acause),]
cvd_cause_ids <- cvd_cause_ids[cause_id %ni% c("ID","ID","ID","ID","ID","ID")]

# Pull data and create three factors for each year (start and end) -------------
loadData <- function() {
  # pull most detailed cause and rei
  cause_ids <- get_cause_metadata(cause_set_id=2, release_id = release_id) %>%
    .[most_detailed==1, cause_id]
  
  ## subset to cvd_cause_ids
  cvd_cause_ids <- cause_ids[cause_ids %in% cvd_cause_ids$cause_id]
  
  risk_draws <- get_draws(gbd_id = cvd_cause_ids, location_id = location_id, year_id = year_ids, sex_id = c(1,2),
                          release_id = release_id, gbd_id_type = 'cause_id', version_id = burdenator_id,
                          measure_id = measure_id, metric_id = c(1,2), source = 'burdenator', age_group_id = age_group_ids,
                          downsample = T, n_draws = ndraws)
  ## parent CVD
  risk_draws2 <- get_draws(gbd_id = "ID", location_id = location_id, year_id = year_ids, sex_id = c(1,2),
                           release_id = release_id, gbd_id_type = 'cause_id', version_id = burdenator_id,
                           measure_id = measure_id, metric_id = c(1,2), source = 'burdenator', age_group_id = age_group_ids,
                           downsample = T, n_draws = ndraws)
  risk_draws <- rbind(risk_draws, risk_draws2)
  
  ## paf
  paf_df <- risk_draws[metric_id == 2,]
  paf_df_long <- melt(paf_df, measure.vars = names(paf_df)[grepl("draw", names(paf_df))], variable.name = "draw", value.name = "met_2")
  paf_df_long <- data.table(paf_df_long)
  
  ## number attr
  num_df <- risk_draws[metric_id == 1,]
  num_df_long <- melt(num_df, measure.vars = names(num_df)[grepl("draw", names(num_df))], variable.name = "draw", value.name = "met_1")
  num_df_long <- data.table(num_df_long)
  
  ## merge together.
  paf_num_df_long <- merge(num_df_long[,c('age_group_id','location_id','sex_id','year_id','cause_id','rei_id','draw','met_1')],
                           paf_df_long[,c('age_group_id','location_id','sex_id','year_id','cause_id','rei_id','draw','met_2')],
                           by=c('age_group_id','location_id','sex_id','year_id','cause_id','rei_id','draw'))
  ## burden
  cdf <- get_draws(gbd_id = cvd_cause_ids, location_id = location_id, year_id = year_ids, sex_id = c(1,2),
                   release_id= release_id, gbd_id_type = 'cause_id', age_group_id = age_group_ids, version_id = dalynator_id,
                   measure_id = measure_id, metric_id=1, source = 'dalynator',
                   downsample = T, n_draws = ndraws)
  
  # for parent CVD
  cdf2 <- get_draws(gbd_id = 'ID', location_id = location_id, year_id = year_ids, sex_id = c(1,2),
                    release_id= release_id, gbd_id_type = 'cause_id', age_group_id = age_group_ids, version_id = dalynator_id,
                    measure_id = measure_id, metric_id=1, source = 'dalynator',
                    downsample = T, n_draws = ndraws)
  
  cdf <- rbind(cdf, cdf2)
  cdf_long <- melt(cdf, measure.vars = names(cdf)[grepl("draw", names(cdf))], variable.name = "draw", value.name = "met_3")
  cdf_long <- data.table(cdf_long)
  
  ## Load pop data
  popdf <- get_population(location_id = location_id, year_id = unique(cdf_long$year_id), age_group_id = unique(cdf_long$age_group_id),
                          sex_id = unique(cdf_long$sex_id), release_id = release_id, run_id = population_estimate_version)
  ## merge on by draw.
  cdf_long <- merge(cdf_long, popdf, by=c('location_id','year_id','sex_id','age_group_id'))
  cdf_long <- data.table(cdf_long)
  cdf_long[,met_3 := met_3/population]
  
  ## Squaring code 
  years <- unique(paf_df_long$year_id)
  sex_id <- unique(paf_df_long$sex_id)
  draws <- paste0('draw_',0:(ndraws-1)) 
  dummy_rei_list <- list()
  for(cause in unique(num_df_long$cause_id)){
    for(rei in unique(num_df_long[cause_id==cause,rei_id])){
      uniq_ages <- unique(num_df_long[cause_id==cause & rei_id==rei, age_group_id])
      sq_ages <- setdiff(c(8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,235),uniq_ages)
      dummy_df <- data.table(expand.grid(age_group_id = sq_ages, year_id = years, sex_id =  sex_id,draw = draws))
      dummy_df[,`:=` (location_id = location_id, cause_id = cause, rei_id = rei, met_1 = 0, met_2 = 0)]
      cause_rei <- paste0(cause,'_',rei)
      dummy_rei_list[[cause_rei]] <- dummy_df
    }
  }
  
  dummy_rei_all <- data.table(rbindlist(dummy_rei_list))
  
  ## append this onto the end
  rfdf <- data.table(rbind(paf_num_df_long,dummy_rei_all))

  # Pull SEVs for SBP, FPG, IKF PAFs of 1
  map <- fread("PATH/cfr_sev.csv")
  
  # Make requisite transformations
  id_cols <- c("location_id", "year_id", "age_group_id", "sex_id", "cause_id")
  
  ## merge
  rfdf <- merge(rfdf, cdf_long[, c(id_cols, "met_3",'draw'), with = FALSE], by = c(id_cols,'draw'), all.x=TRUE)
  
  ## remove the individual components, no longer needed.
  rm(paf_df, num_df, cdf, sevdf, paf_df_long, num_df_long, cdf_long)
  
  ## Fill na's with 0
  rfdf[is.na(met_1), met_1 := 0] ## met_1: attributable deaths
  rfdf[is.na(met_2), met_2 := 0] ## met_2: PAFs
  rfdf[is.na(met_3), met_3 := 0] ## met_3: CSMR for the cause.
  
  # Flag PAFs of 1 as needed
  paf1_strats <- fread("PATH/strategies_for_pafsofone.csv") %>%  
    .[, .(rei_id, cause_id, paf1)]
  paf1_strats <- merge(paf1_strats,reis[, .(rei_id,parent_id)],by="rei_id")
  paf1_strats[, rei_id := NULL]
  setnames(paf1_strats, "parent_id","rei_id")
  paf1_strats <- unique(paf1_strats)
  rfdf <- merge(rfdf, paf1_strats, by = c("rei_id", "cause_id"), all.x=T)
  rfdf[is.na(paf1), paf1 := "no"]
  
  # Calculate burden factors (exposure and risk-deleted)
  rfdf[, `:=` (under_rate=met_3 * (1 - met_2), ## risk deleted = csmr * (1-paf)
               risk=met_2 / (1 - met_2),       ## risk exposure = paf/(1-paf)
               number=met_1)]                  ## # of attributable deaths
  rfdf[paf1 == "two_risk" | paf1 == "two_under", risk := met_3 * met_2]
  rfdf[paf1 == "cfr_sev" | paf1 == "cfr_prev", under_rate := met_3 ]
  rfdf <- rfdf[!is.na(under_rate) & !is.na(risk)]
  
  # backup testing
  rfdf_backup <- copy(rfdf)
  
  # retrieve population
  popdf <- get_population(location_id = location_id, year_id = year_ids, age_group_id = age_group_ids, 
                            sex_id = sex_ids, release_id = release_id, run_id = population_estimate_version)
  # Combine to one data frame
  id_cols <- id_cols[id_cols != 'cause_id'] 
  df <- merge(
    rfdf[, c(id_cols, "cause_id", "rei_id","draw", "paf1", "under_rate", "risk", "number"), with = FALSE], 
    popdf[, c(id_cols, "population"), with = FALSE], by = c(id_cols))
  
  ## calculate in cases of SEV or PREV
  df[paf1 == "cfr_sev", `:=` (risk=sev,
                              under_rate=under_rate/sev)]
  df[paf1 == "cfr_prev", `:=` (risk=prev/(population),
                               under_rate=under_rate/(prev/(population)))]
  df[is.na(risk), risk := 0]
  df[is.na(under_rate), under_rate := 0]
  df[paf1 == "two_under" & number == 0, under_rate := 0]
  
  # Make wide
  df[, year_id := as.numeric(factor(year_id))]
  df <- data.table::dcast(df, location_id + age_group_id + sex_id + rei_id + cause_id + paf1 + draw ~ year_id,
                          value.var = c("population", "under_rate", "risk", "number"))
  # Return data frame
  return(df)
}

# Run das gupta to calc three effects ------------------------------------------
decompFactors <- function(df_) {
  df <- copy(df_)
  
  # Calculate 3-factor effects
  df[, population_effect := ((under_rate_1 * risk_1 + under_rate_2 * risk_2) / 3 +
                               (under_rate_1 * risk_2 + under_rate_2 * risk_1) / 6) * (population_2 - population_1)]
  df[, under_rate_effect := ((population_1 * risk_1 + population_2 * risk_2) / 3 +
                               (population_1 * risk_2 + population_2 * risk_1) / 6) * (under_rate_2 - under_rate_1)]
  df[, risk_effect := ((population_1 * under_rate_1 + population_2 * under_rate_2) / 3 +
                         (population_1 * under_rate_2 + population_2 * under_rate_1) / 6) * (risk_2 - risk_1)]
  
  # Calculate 2-factor effects for PAFs of 1
  df[paf1 == "two_risk" | paf1 == "two_under", `:=` (
    population_effect=((risk_1 + risk_2) / 2) * (population_2 - population_1),
    risk_effect=((population_1 + population_2) / 2) * (risk_2 - risk_1))]
  
  # set to 0 as needed
  df[paf1 == "two_risk", under_rate_effect := 0]
  df[paf1 == "two_under", under_rate_effect := risk_effect]
  df[paf1 == "two_under", risk_effect := 0]
  df[is.na(risk_effect), risk_effect := 0]
  df[is.na(population_effect), population_effect := 0]
  df[is.na(under_rate_effect), under_rate_effect := 0]
  df[number_1==0 & number_2==0, `:=` (risk_effect=0,
                                      under_rate_effect=0,
                                      population_effect=0)]
  return(df)
}
# Apply mediation --------------------------------------------------------------
meatiator <- function(df, reimetadf) {
  # Get mediation data
  
  meddf <- fread("PATH/mediation_matrix_draw_gbd_2023.csv")
  meddf[, mean := rowMeans(.SD, na.rm=T), .SDcols=paste0("draw_", 0:999)]
  meddf <- meddf[, grep("draw", names(meddf), value = TRUE) := .(rep(NULL, 1000))][!is.na(mean),]
  
  # Mediate change factor
  overlapdf <- merge(df,
                     meddf[, c("rei_id", "cause_id", "mean", "med_id"), with = FALSE],
                     by = c("rei_id", "cause_id"))
  overlapdf[, risk_effect_mediation := risk_effect * mean]
  
  overlapdf <- merge(overlapdf[, c("location_id", "cause_id", "sex_id", "age_group_id", "draw",
                                   "med_id", "risk_effect_mediation"), with = FALSE],
                     reimetadf[, c("rei_id", "rei"), with = FALSE],
                     by.x = "med_id",
                     by.y = "rei_id")
  overlapdf[,rei_id := med_id]
  overlapdf <- overlapdf[, .(risk_effect_mediation = sum(risk_effect_mediation)),
                         by = .(location_id, rei_id, cause_id, sex_id, age_group_id,draw)]
  # Apply mediation
  df <- merge(df, overlapdf,
              by = c("location_id", "cause_id", "rei_id", "sex_id","age_group_id",'draw'),
              all.x = TRUE)
  
  ## set to 0 in case of contrasting directions
  df[(risk_effect < 0 & risk_effect_mediation > 0) | (risk_effect > 0 & risk_effect_mediation < 0), risk_effect_mediation := 0] 
  df[abs(risk_effect) < abs(risk_effect_mediation), risk_effect_mediation := risk_effect]
  df[!is.na(risk_effect_mediation), risk_effect := risk_effect - risk_effect_mediation]
  df[, risk_effect_mediation := NULL]
  rm(meddf, overlapdf)
  # Return data frame
  return(df)
}
# Scale mediated most detailed risks to all-risk aggregate ---------------------
riskRaker <- function(df, reimetadf) {
  ## test for the under-effect ##
  df[, agg_0 := ((number_2-number_1)-(risk_effect+population_effect))]
  df[under_rate_effect >= 0, `:=` (agg_1=under_rate_effect, agg_2=0)]
  df[under_rate_effect < 0, `:=` (agg_1=0, agg_2=under_rate_effect)]
  df[,scalar_1:=(agg_0+sqrt(agg_0^2-4*agg_1*agg_2))/(2*agg_1)]
  df[,scalar_2:=1/((agg_0+sqrt(agg_0^2-4*agg_1*agg_2))/(2*agg_1))]
  df[agg_0==0 & (agg_1==0|agg_2==0), `:=` (scalar_1=0, scalar_2=0)]
  df[agg_0<0 & agg_1==0 & agg_2!=0, `:=` (scalar_1=1, scalar_2=agg_0/agg_2)]
  df[agg_0>0 & agg_2==0 & agg_1!=0, `:=` (scalar_1=agg_0/agg_1, scalar_2=1)]
  df[agg_0>0 & agg_1==0 & agg_2!=0, `:=` (scalar_1=1, scalar_2=agg_0/agg_2)]
  df[agg_0<0 & agg_2==0 & agg_1!=0, `:=` (scalar_1=agg_0/agg_1, scalar_2=1)]
  df[agg_0!=0 & agg_1==0 & agg_2==0, `:=` (scalar_1=0, scalar_2=0)]
  
  # make long to prep for merge back on scalars for negative and positive values
  df[under_rate_effect <0, under_rate_effect := under_rate_effect * scalar_2]
  df[under_rate_effect >= 0, under_rate_effect := under_rate_effect * scalar_1]
  df[, c("agg_1","agg_2", "agg_0", "scalar_1", "scalar_2","factor_scalar") := .(NULL,NULL, NULL, NULL, NULL,NULL)]
  # Return data frame
  return(df)
}
# Run! -------------------------------------------------------------------------
# pull data and run decomp
message('DONE ESTABLISHING FUNCTIONS')
df_data_ <- loadData()
message('DONE GETTING DATA')
df_factor <- decompFactors(df_data_)
message('DONE COMPUTING FACTORS')
df_med <- meatiator(df_factor, reimetadf)
message('DONE MEDIATING')
df_rake <- riskRaker(df_med, reimetadf)
message('DONE RAKING')

## aggregate
df_final_med <- df_rake[, lapply(.SD, sum), by = c("location_id", "rei_id", "cause_id","draw"),
                        .SDcols = c("number_1", "number_2", "population_effect", "under_rate_effect", "risk_effect")]
df_final_med <- df_final_med[rei_id %in% reimetadf$rei_id]

measure <- ifelse(measure_id == 2, "daly", ifelse(measure_id == 1, "death", ifelse(measure_id==2, "yld", "yll")))
meas_name <- ifelse(measure_id == 2, "DALY",
                    ifelse(measure_id == 1, "mortality",
                           ifelse(measure_id==2, "YLD", "YLL")))
meas_names <- ifelse(measure_id == 1, meas_name, paste0(meas_name, "s"))


# calculate percents
loc_set <- copy(location_id)
df_final_med <- df_final_med[location_id == loc_set,][,location_id := NULL]
df_final_med[, sex_id := 3]
df_final_med[, total_pct := (number_2 - number_1) / number_1]

## save number change for later.
df_number_change <- copy(df_final_med)
df_number_change[,total_burden_change := (number_2-number_1)]

## numbers get reassigned to percents
df_final_med[, c("population_effect_pct", "under_rate_effect_pct", "risk_effect_pct") :=
               .(population_effect / number_1,
                 under_rate_effect / number_1,
                 risk_effect / number_1)]

## save counts for later
df_final_med_count <- copy(df_final_med[,c('sex_id','rei_id','cause_id','draw','number_1','number_2','population_effect','under_rate_effect','risk_effect','total_pct')])

## remove unnecessary columns and rename
df_final_med[, c("number_1", "number_2", "population_effect","under_rate_effect","risk_effect") := .(NULL, NULL, NULL, NULL, NULL)]
setnames(df_final_med, c("population_effect_pct","under_rate_effect_pct","risk_effect_pct"),c("population_effect","under_rate_effect","risk_effect"))

# retrieve population draws for uncertainty in change due to aging and growth
pop_draws_path <- paste0("PATH",population_estimate_version,"/",location_id,"PATH/population_reporting_raked_ui.h5")
pop_list <- list()
for(draw in 0:(ndraws-1)){
  poptest <- data.table(h5read(pop_draws_path, name=as.character(draw)))
  poptest <- poptest[age_group_id %in% 22 & year_id %in% year_ids & sex_id==3,]
  pop_list[[as.character(draw)]] <- poptest
}
pop <- data.table(rbindlist(pop_list))
pop[,draw := paste0('draw_',draw)]
setnames(pop, c('value'),c('population'))
pop <- as.data.table(dcast(pop,location_id + sex_id + draw ~ year_id, value.var="population"))
setnames(pop, paste0(year_ids), paste0("yr", year_ids))
pop[, popchg := (get(paste0("yr",year_ids[2])) - get(paste0("yr",year_ids[1]))) / get(paste0("yr",year_ids[1]))]

## for percent
df_final_med <- merge(df_final_med, pop[,.(sex_id,popchg,draw)], by=c("sex_id",'draw'))
df_final_med[, age_structure_effect := population_effect]
df_final_med[, population_effect := popchg]

# age structure effect
df_final_med[, age_structure_effect := age_structure_effect - popchg]
df_final_med[, popchg:= NULL]
df_final_med <- data.table::melt(df_final_med, id.vars = c("cause_id", "rei_id", "sex_id","draw"),
                                 value.vars = c(grep("effect", names(df_final_med), value = TRUE), "total_pct"),
                                 variable.name = "factor_format", value.name = "change")
## for count
df_final_med_count <- merge(df_final_med_count, pop[,.(sex_id,popchg,draw)], by=c("sex_id","draw"))

## directly applying percentages
df_final_med_count <- merge(df_final_med_count,
                            df_final_med[factor_format=='age_structure_effect',c('cause_id','rei_id','sex_id','draw','change')],
                            by = c('cause_id','rei_id','sex_id','draw'))
setnames(df_final_med_count, 'change','change_age')
df_final_med_count <- merge(df_final_med_count,
                            df_final_med[factor_format=='population_effect',c('cause_id','rei_id','sex_id','draw','change')],
                            by = c('cause_id','rei_id','sex_id','draw'))
setnames(df_final_med_count, 'change','change_pop')
df_final_med_count[, population_effect := number_1*change_pop]
df_final_med_count[, age_structure_effect := number_1*change_age]
df_final_med_count[, `:=` (change_pop = NULL, change_age = NULL)]
df_final_med_count <- df_final_med_count[!(is.na(total_pct) | is.infinite(total_pct) | is.nan(total_pct)),]

## Calculate percents from the count.
df_final_med_mean <- df_final_med_count[,.(mean_n1 = mean(number_1),
                                           mean_n2 = mean(number_2),
                                           mean_pop_e = mean(population_effect),
                                           mean_und_e = mean(under_rate_effect),
                                           mean_ris_e = mean(risk_effect),
                                           mean_age_e = mean(age_structure_effect)),
                                        by = c('rei_id','cause_id','sex_id')]

## calculate percent at the mean level.
df_final_med_mean <- df_final_med_mean[,c("population_effect_pct", "under_rate_effect_pct", "risk_effect_pct",'age_structure_effect_pct','total_pct') :=
                                         .(mean_pop_e / mean_n1,
                                           mean_und_e / mean_n1,
                                           mean_ris_e / mean_n1,
                                           mean_age_e / mean_n1,
                                           (mean_n2 - mean_n1) / mean_n1)]

## calculate uncertainty for percent and count
df_final_med_count[,c("population_effect_pct_draw", "under_rate_effect_pct_draw", "risk_effect_pct_draw", "age_structure_effect_pct_draw", 'total_pct_draw') :=
                     .(population_effect / number_1,
                       under_rate_effect / number_1,
                       risk_effect / number_1,
                       age_structure_effect / number_1,
                       (number_2 - number_1) / number_1)]
df_final_med_count <- df_final_med_count[!is.na(population_effect_pct_draw)]

## percent space
df_final_med_pct_uncert <- df_final_med_count[,.(lower_population_effect_pct = quantile(population_effect_pct_draw, 0.025),
                                                 upper_population_effect_pct = quantile(population_effect_pct_draw, 0.975),
                                                 lower_under_rate_effect_pct = quantile(under_rate_effect_pct_draw, 0.025),
                                                 upper_under_rate_effect_pct = quantile(under_rate_effect_pct_draw, 0.975),
                                                 lower_risk_effect_pct = quantile(risk_effect_pct_draw, 0.025),
                                                 upper_risk_effect_pct = quantile(risk_effect_pct_draw, 0.975),
                                                 lower_age_structure_effect_pct = quantile(age_structure_effect_pct_draw, 0.025),
                                                 upper_age_structure_effect_pct = quantile(age_structure_effect_pct_draw, 0.975),
                                                 lower_total_pct = quantile(total_pct_draw, 0.025),
                                                 upper_total_pct = quantile(total_pct_draw, 0.975)),
                                              by = c('rei_id','cause_id','sex_id')]

## number space
df_final_med_count_uncert <- df_final_med_count[,.(lower_population_effect = quantile(population_effect, 0.025),
                                                   upper_population_effect = quantile(population_effect, 0.975),
                                                   lower_under_rate_effect = quantile(under_rate_effect, 0.025),
                                                   upper_under_rate_effect = quantile(under_rate_effect, 0.975),
                                                   lower_risk_effect = quantile(risk_effect, 0.025),
                                                   upper_risk_effect = quantile(risk_effect, 0.975),
                                                   lower_age_structure_effect = quantile(age_structure_effect, 0.025),
                                                   upper_age_structure_effect = quantile(age_structure_effect, 0.975)),
                                                by = c('rei_id','cause_id','sex_id')]
# reformat data - mean percent #
df_final_med <- df_final_med_mean[,c('rei_id','cause_id','sex_id','population_effect_pct','age_structure_effect_pct','under_rate_effect_pct','risk_effect_pct','total_pct')]

setnames(df_final_med, c("population_effect_pct","under_rate_effect_pct","risk_effect_pct",'age_structure_effect_pct'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect'))

df_final_med <- data.table::melt(df_final_med, id.vars = c("cause_id", "rei_id", "sex_id"),
                                 value.vars = c(grep("effect", names(df_final_med), value = TRUE), "total_pct"),
                                 variable.name = "factor_format", value.name = "change")
total_df_final_med <- df_final_med[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med, "change", "total")

df_final_med <- merge(df_final_med[factor_format!="total_pct",], total_df_final_med, by=c("cause_id","rei_id","sex_id"))
df_final_med <- merge(df_final_med, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level.
risk_level <- 2
df_final_med <- df_final_med[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 
df_final_med[, factor_format := factor(factor_format,
                                       level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                 "risk_effect"),
                                       label = c("Change due to population growth",
                                                 "Change due to population ageing",
                                                 paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                 "Change due to risk exposure"))]
df_final_med[, location_id := loc_set]
df_final_med[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## calculation upper/lower bounds for percents ## 

## subset to CVD causes - pct lower bound calculations
df_final_med <- df_final_med[cause_id %in% cvd_cause_ids$cause_id]
df_final_med_ui_lower <- df_final_med_pct_uncert[,c('rei_id','cause_id','sex_id',
                                                    'lower_population_effect_pct',
                                                    'lower_age_structure_effect_pct',
                                                    'lower_under_rate_effect_pct',
                                                    'lower_risk_effect_pct',
                                                    'lower_total_pct')]
setnames(df_final_med_ui_lower, c("lower_population_effect_pct","lower_under_rate_effect_pct","lower_risk_effect_pct",'lower_age_structure_effect_pct','lower_total_pct'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect','total_pct'))
df_final_med_ui_lower <- data.table::melt(df_final_med_ui_lower, id.vars = c("cause_id", "rei_id", "sex_id"),
                                          value.vars = c(grep("effect", names(df_final_med_ui_lower), value = TRUE), "total_pct"),
                                          variable.name = "factor_format", value.name = "change")
total_df_final_med_ui_lower <- df_final_med_ui_lower[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med_ui_lower, "change", "total")
df_final_med_ui_lower <- merge(df_final_med_ui_lower[factor_format!="total_pct",], total_df_final_med_ui_lower, by=c("cause_id","rei_id","sex_id"))
df_final_med_ui_lower <- merge(df_final_med_ui_lower, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level - pct lower bound calculations
risk_level <- 2
df_final_med_ui_lower <- df_final_med_ui_lower[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 
df_final_med_ui_lower[, factor_format := factor(factor_format,
                                                level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                          "risk_effect"),
                                                label = c("Change due to population growth",
                                                          "Change due to population ageing",
                                                          paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                          "Change due to risk exposure"))]
df_final_med_ui_lower[, location_id := loc_set]
df_final_med_ui_lower[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## subset to CVD causes - pct upper bound calculations
df_final_med_ui_lower <- df_final_med_ui_lower[cause_id %in% cvd_cause_ids$cause_id]
setnames(df_final_med_ui_lower, c('change','total'), c('lower_change','lower_total'))
df_final_med_ui_upper <- df_final_med_pct_uncert[,c('rei_id','cause_id','sex_id',
                                                    'upper_population_effect_pct',
                                                    'upper_age_structure_effect_pct',
                                                    'upper_under_rate_effect_pct',
                                                    'upper_risk_effect_pct',
                                                    'upper_total_pct')]
setnames(df_final_med_ui_upper, c("upper_population_effect_pct","upper_under_rate_effect_pct","upper_risk_effect_pct",'upper_age_structure_effect_pct','upper_total_pct'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect','total_pct'))
df_final_med_ui_upper <- data.table::melt(df_final_med_ui_upper, id.vars = c("cause_id", "rei_id", "sex_id"),
                                          value.vars = c(grep("effect", names(df_final_med_ui_upper), value = TRUE), "total_pct"),
                                          variable.name = "factor_format", value.name = "change")
total_df_final_med_ui_upper <- df_final_med_ui_upper[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med_ui_upper, "change", "total")
df_final_med_ui_upper <- merge(df_final_med_ui_upper[factor_format!="total_pct",], total_df_final_med_ui_upper, by=c("cause_id","rei_id","sex_id"))
df_final_med_ui_upper <- merge(df_final_med_ui_upper, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level - pct upper calculations
risk_level <- 2
df_final_med_ui_upper <- df_final_med_ui_upper[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 
df_final_med_ui_upper[, factor_format := factor(factor_format,
                                                level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                          "risk_effect"),
                                                label = c("Change due to population growth",
                                                          "Change due to population ageing",
                                                          paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                          "Change due to risk exposure"))]
df_final_med_ui_upper[, location_id := loc_set]
df_final_med_ui_upper[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## subset to CVD causes - number mean
df_final_med_ui_upper <- df_final_med_ui_upper[cause_id %in% cvd_cause_ids$cause_id]
setnames(df_final_med_ui_upper, c('change','total'), c('upper_change','upper_total'))

# reformat data - count #
df_final_med_count_mean <- df_final_med_mean[,c('rei_id','cause_id','sex_id','mean_pop_e','mean_und_e','mean_ris_e','mean_age_e')]

setnames(df_final_med_count_mean, c('mean_pop_e','mean_und_e','mean_ris_e','mean_age_e'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect'))
df_final_med_count_mean <- data.table::melt(df_final_med_count_mean, id.vars = c("cause_id", "rei_id", "sex_id"),
                                            value.vars = c(grep("effect", names(df_final_med_count_mean), value = TRUE), "total_pct"),
                                            variable.name = "factor_format", value.name = "change")
total_df_final_med_count_mean <- df_final_med_count_mean[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med_count_mean, "change", "total")
df_final_med_count_mean <- merge(df_final_med_count_mean[factor_format!="total_pct",], total_df_final_med_count_mean, by=c("cause_id","rei_id","sex_id"))
df_final_med_count_mean <- merge(df_final_med_count_mean, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level  - number mean
risk_level <- 2
df_final_med_count_mean <- df_final_med_count_mean[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 
df_final_med_count_mean[, factor_format := factor(factor_format,
                                                  level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                            "risk_effect"),
                                                  label = c("Change due to population growth",
                                                            "Change due to population ageing",
                                                            paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                            "Change due to risk exposure"))]
df_final_med_count_mean[, location_id := loc_set]
df_final_med_count_mean[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## Calculate lower/upper bound for counts ## 
## subset to CVD causes - 
df_final_med_count_mean <- df_final_med_count_mean[cause_id %in% cvd_cause_ids$cause_id]

# reformat data - count #
df_final_med_count_mean <- df_final_med_mean[,c('rei_id','cause_id','sex_id','mean_pop_e','mean_und_e','mean_ris_e','mean_age_e')]
setnames(df_final_med_count_mean, c('mean_pop_e','mean_und_e','mean_ris_e','mean_age_e'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect'))
df_final_med_count_mean <- data.table::melt(df_final_med_count_mean, id.vars = c("cause_id", "rei_id", "sex_id"),
                                            value.vars = c(grep("effect", names(df_final_med_count_mean), value = TRUE), "total_pct"),
                                            variable.name = "factor_format", value.name = "change")
total_df_final_med_count_mean <- df_final_med_count_mean[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]

# reformat data - count #
df_final_med_count_mean <- df_final_med_mean[,c('rei_id','cause_id','sex_id','mean_pop_e','mean_und_e','mean_ris_e','mean_age_e','total_pct')]
setnames(df_final_med_count_mean, c('mean_pop_e','mean_und_e','mean_ris_e','mean_age_e'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect'))
df_final_med_count_mean <- data.table::melt(df_final_med_count_mean, id.vars = c("cause_id", "rei_id", "sex_id"),
                                            value.vars = c(grep("effect", names(df_final_med_count_mean), value = TRUE), "total_pct"),
                                            variable.name = "factor_format", value.name = "change")
total_df_final_med_count_mean <- df_final_med_count_mean[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med_count_mean, "change", "total")
df_final_med_count_mean <- merge(df_final_med_count_mean[factor_format!="total_pct",], total_df_final_med_count_mean, by=c("cause_id","rei_id","sex_id"))
df_final_med_count_mean <- merge(df_final_med_count_mean, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level.
risk_level <- 2
df_final_med_count_mean <- df_final_med_count_mean[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 

df_final_med_count_mean[, factor_format := factor(factor_format,
                                                  level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                            "risk_effect"),
                                                  label = c("Change due to population growth",
                                                            "Change due to population ageing",
                                                            paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                            "Change due to risk exposure"))]
df_final_med_count_mean[, location_id := loc_set]
df_final_med_count_mean[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## subset to CVD causes
df_final_med_count_mean <- df_final_med_count_mean[cause_id %in% cvd_cause_ids$cause_id]

df_final_med_count_lower <- df_final_med_count_uncert[,c('rei_id','cause_id','sex_id','lower_population_effect','lower_under_rate_effect','lower_risk_effect','lower_age_structure_effect')]
df_final_med_count_lower[,total_pct := NA]
df_final_med_count_lower[,total_pct := NA] 
setnames(df_final_med_count_lower, c('lower_population_effect','lower_under_rate_effect','lower_risk_effect','lower_age_structure_effect'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect'))
df_final_med_count_lower <- data.table::melt(df_final_med_count_lower, id.vars = c("cause_id", "rei_id", "sex_id"),
                                             value.vars = c(grep("effect", names(df_final_med_count_lower), value = TRUE), "total_pct"),
                                             variable.name = "factor_format", value.name = "change")
total_df_final_med_count_lower <- df_final_med_count_lower[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med_count_lower, "change", "total")
df_final_med_count_lower <- merge(df_final_med_count_lower[factor_format!="total_pct",], total_df_final_med_count_lower, by=c("cause_id","rei_id","sex_id"))
df_final_med_count_lower <- merge(df_final_med_count_lower, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level.
risk_level <- 2
df_final_med_count_lower <- df_final_med_count_lower[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 
df_final_med_count_lower[, factor_format := factor(factor_format,
                                                   level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                             "risk_effect"),
                                                   label = c("Change due to population growth",
                                                             "Change due to population ageing",
                                                             paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                             "Change due to risk exposure"))]
df_final_med_count_lower[, location_id := loc_set]
df_final_med_count_lower[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## subset to CVD causes
df_final_med_count_lower <- df_final_med_count_lower[cause_id %in% cvd_cause_ids$cause_id]

# reformat data - upper UI lower#
df_final_med_count_upper <- df_final_med_count_uncert[,c('rei_id','cause_id','sex_id','upper_population_effect','upper_under_rate_effect','upper_risk_effect','upper_age_structure_effect')]
df_final_med_count_upper[,total_pct := NA] 
setnames(df_final_med_count_upper, c('upper_population_effect','upper_under_rate_effect','upper_risk_effect','upper_age_structure_effect'),
         c("population_effect","under_rate_effect","risk_effect",'age_structure_effect'))
df_final_med_count_upper <- data.table::melt(df_final_med_count_upper, id.vars = c("cause_id", "rei_id", "sex_id"),
                                             value.vars = c(grep("effect", names(df_final_med_count_upper), value = TRUE), "total_pct"),
                                             variable.name = "factor_format", value.name = "change")
total_df_final_med_count_upper <- df_final_med_count_upper[factor_format=="total_pct", list(cause_id,rei_id,change,sex_id)]
setnames(total_df_final_med_count_upper, "change", "total")
df_final_med_count_upper <- merge(df_final_med_count_upper[factor_format!="total_pct",], total_df_final_med_count_upper, by=c("cause_id","rei_id","sex_id"))
df_final_med_count_upper <- merge(df_final_med_count_upper, reimetadf[, c("rei_id", "rei_name", "sort_order", "level", "most_detailed"), with = FALSE], by = "rei_id")

## set risk level.
risk_level <- 2
df_final_med_count_upper <- df_final_med_count_upper[level %in% c(0,1,risk_level) | (most_detailed == 1 & level < risk_level)] 
df_final_med_count_upper[, factor_format := factor(factor_format,
                                                   level = c("population_effect", "age_structure_effect", "under_rate_effect",
                                                             "risk_effect"),
                                                   label = c("Change due to population growth",
                                                             "Change due to population ageing",
                                                             paste0("Change due to risk-deleted ", meas_name, " rate"),
                                                             "Change due to risk exposure"))]
df_final_med_count_upper[, location_id := loc_set]
df_final_med_count_upper[, location_name := unique(get_location_metadata(location_set_id=35, release_id = release_id
)[level > 3,location_name := paste0(
  lancet_label," (", substr(ihme_loc_id, 1, 3),")")] %>%
  .[location_id==loc_set,] %>% .$location_name)]

## subset to CVD causes
df_final_med_count_upper <- df_final_med_count_upper[cause_id %in% cvd_cause_ids$cause_id]


df_final_med_pct <- merge(df_final_med,
                          df_final_med_ui_lower[,c('rei_id','location_id','cause_id','sex_id','factor_format','lower_change','lower_total')],
                          by=c('rei_id','location_id','cause_id','sex_id','factor_format'))

df_final_med_pct <- merge(df_final_med_pct,
                          df_final_med_ui_upper[,c('rei_id','location_id','cause_id','sex_id','factor_format','upper_change','upper_total')],
                          by=c('rei_id','location_id','cause_id','sex_id','factor_format'))


setnames(df_final_med_count_mean, c('change','total'), c('change_count','count_total'))
setnames(df_final_med_count_lower, c('change','total'), c('lower_change_count','lower_count_total'))
setnames(df_final_med_count_upper, c('change','total'), c('upper_change_count','upper_count_total'))

## merge together count lower/upper
df_final_med_count_withui <- merge(df_final_med_count_mean,
                                   df_final_med_count_lower[,c('rei_id','location_id','cause_id','sex_id','factor_format','lower_change_count')],
                                   by = c('rei_id','location_id','cause_id','sex_id','factor_format'))
df_final_med_count_withui <- merge(df_final_med_count_withui,
                                   df_final_med_count_upper[,c('rei_id','location_id','cause_id','sex_id','factor_format','upper_change_count')],
                                   by = c('rei_id','location_id','cause_id','sex_id','factor_format'))


## merge pct and count
df_final_med_save <- merge(df_final_med_pct[,c('rei_id','location_id','cause_id','sex_id','factor_format','change','lower_change','upper_change','total','lower_total','upper_total')],
                           df_final_med_count_withui[,c('rei_id','location_id','cause_id','sex_id','factor_format','change_count','lower_change_count','upper_change_count')],
                           by = c('rei_id','location_id','cause_id','sex_id','factor_format'))


openxlsx::write.xlsx(df_final_med_save,paste0(outdir,'measure_',measure_id,'/',location_id,'_decomp.xlsx'))
